"""Read, analyze, and visualize molecular structures from XYZ files.

It includes a `MolGraph` class that represents a molecular graph structure, with methods for
setting element radii and colors, filtering atoms, reading XYZ files, generating adjacency lists,
and converting to Plotly and NetworkX representations.
"""

import re
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import (
    Dict,
    FrozenSet,
    Iterator,
    List,
    Optional,
    Sequence,
    Set,
    Tuple,
    Union,
)

import networkx as nx
import numpy as np
import plotly.graph_objs as go
from numpy.typing import NDArray

from xyz2graph.geometry import Point3D
from xyz2graph.molecule import Atom, Bond
from xyz2graph.visualization import VisualizationConfig, create_visualization

from .constants import _DEFAULT_CPK_COLORS, _DEFAULT_RADII
from .logging import logger


@dataclass
class MolGraph:
    """Represents a molecular graph structure.

    This class provides functionality to read molecular structure data from XYZ files,
    analyze molecular geometry, and visualize the molecular structure using Plotly.

    Examples:
        >>> mg = MolGraph()
        >>> mg.read_xyz('molecule.xyz')
        >>> print(mg.formula())
        'C2H6O'
    """

    atoms: List[Atom] = field(default_factory=list)
    bonds: List[Bond] = field(default_factory=list)
    comment: str = field(default="")

    # Customizable parameters with defaults
    default_radii: Dict[str, float] = field(default_factory=lambda: dict(_DEFAULT_RADII))
    cpk_colors: Dict[str, str] = field(default_factory=lambda: dict(_DEFAULT_CPK_COLORS))
    cpk_color_rest: str = field(default="pink")

    @property
    def elements(self) -> List[str]:
        """Get list of elements in the molecule."""
        return [atom.element for atom in self.atoms]

    @property
    def x(self) -> List[float]:
        """Get list of x coordinates."""
        return [atom.x for atom in self.atoms]

    @property
    def y(self) -> List[float]:
        """Get list of y coordinates."""
        return [atom.y for atom in self.atoms]

    @property
    def z(self) -> List[float]:
        """Get list of z coordinates."""
        return [atom.z for atom in self.atoms]

    @property
    def xyz(self) -> NDArray[np.float64]:
        """Get coordinates as a Nx3 numpy array."""
        return np.array([[atom.x, atom.y, atom.z] for atom in self.atoms])

    @property
    def atomic_radii(self) -> List[float]:
        """Get list of atomic radii."""
        return [atom.radius for atom in self.atoms]

    @property
    def bond_lengths(self) -> Dict[FrozenSet[int], float]:
        """Get dictionary of bond lengths indexed by atom pairs."""
        return {bond.indices: bond.length for bond in self.bonds}

    @property
    def indices(self) -> List[int]:
        """Get list of atom indices."""
        return [atom.index for atom in self.atoms]

    @property
    def adj_list(self) -> Dict[int, Set[int]]:
        """Get adjacency list representation of molecular graph."""
        adj: Dict[int, Set[int]] = {}
        for bond in self.bonds:
            i1, i2 = tuple(bond.indices)
            adj.setdefault(i1, set()).add(i2)
            adj.setdefault(i2, set()).add(i1)
        return adj

    @property
    def adj_matrix(self) -> NDArray[np.int_]:
        """Get adjacency matrix representation of molecular graph.

        Returns:
            NDArray[np.int_]: A square matrix where entry (i,j) is 1 if atoms i and j
            are bonded and 0 otherwise.
        """
        n = len(self.atoms)
        matrix = np.zeros((n, n), dtype=np.int_)
        for bond in self.bonds:
            i, j = tuple(bond.indices)
            matrix[i, j] = matrix[j, i] = 1
        return matrix

    @property
    def distance_matrix(self) -> NDArray[np.float64]:
        """Get matrix of interatomic distances.

        Returns:
            NDArray[np.float64]: A square matrix where entry (i,j)
                is the distance between atoms i and j.
        """
        n_atoms = len(self.atoms)
        MEMORY_THRESHOLD = 10_000  # Switch to loop method above 10k atoms

        try:
            if n_atoms < MEMORY_THRESHOLD:
                # Fast broadcasting method for small/medium molecules
                distances = self.xyz[:, np.newaxis, :] - self.xyz
                return np.sqrt(np.einsum("ijk,ijk->ij", distances, distances))
            else:
                logger.info(
                    f"Large molecule ({n_atoms} atoms) detected. Using memory-efficient method."
                )
        except MemoryError:
            logger.info("Memory allocation failed. Switching to memory-efficient method.")

        # Memory efficient loop method for large molecules or if broadcasting fails
        distance_matrix = np.zeros((n_atoms, n_atoms), dtype=np.float64)
        for i in range(n_atoms):
            diff = self.xyz[i] - self.xyz
            distance_matrix[i] = np.sqrt(np.sum(diff * diff, axis=1))

        return distance_matrix

    def _generate_bonds(self) -> None:
        """Generate bonds based on atomic distances."""
        self.bonds.clear()

        distances = self.distance_matrix

        # Calculate bonding threshold matrix
        radii = np.array(self.atomic_radii)
        thresholds = (radii[:, np.newaxis] + radii) * 1.3

        # Find bonded atoms (creates adjacency matrix internally)
        adj_matrix = np.logical_and(0.1 < distances, thresholds > distances).astype(np.int_)

        # Create Bond objects from adjacency matrix
        for i, j in zip(*np.nonzero(adj_matrix)):
            if i < j:  # Avoid duplicate bonds
                self.bonds.append(Bond(self.atoms[i], self.atoms[j]))

    def set_element_radius(self, element: str, radius: float) -> None:
        """Set the reference radius for a specific element and update bonds."""
        self.default_radii[element] = radius
        # Update radii for existing atoms of this element
        for atom in self.atoms:
            if atom.element == element:
                atom.radius = radius
        # Regenerate bonds with new radii
        self._generate_bonds()

    def filter(
        self,
        indices: Optional[List[int]] = None,
        elements: Optional[List[str]] = None,
        inplace: bool = False,
    ) -> Optional["MolGraph"]:
        """Filter atoms by indices and/or elements."""
        # Create mask for atoms to keep
        mask = [True] * len(self.atoms)

        if indices is not None:
            # Validate indices
            if any(i < 0 or i >= len(self.atoms) for i in indices):
                raise IndexError("Atom index out of range")
            for idx in indices:
                mask[idx] = False

        if elements is not None:
            # Track which filter elements are actually used
            found_elements = set()
            for atom in self.atoms:
                if atom.element in elements:
                    found_elements.add(atom.element)
                    mask[atom.index] = False

            # Warn about unused filter elements
            unused_elements = set(elements) - found_elements
            if unused_elements:
                logger.warning(
                    f"Element(s) not found: {', '.join(sorted(unused_elements))}. "
                    "Use proper case (e.g., 'H' not 'h')"
                )

        if not any(mask):
            raise ValueError("Cannot filter out all atoms from molecule")

        filtered_atoms = [atom for atom, keep in zip(self.atoms, mask) if keep]

        if inplace:
            self.atoms = filtered_atoms
            self._generate_bonds()
            return None

        # Create new instance
        new_mol = MolGraph()
        new_mol.atoms = filtered_atoms
        new_mol.comment = self.comment
        new_mol.default_radii = self.default_radii.copy()
        new_mol.cpk_colors = self.cpk_colors.copy()
        new_mol.cpk_color_rest = self.cpk_color_rest
        new_mol._generate_bonds()
        return new_mol

    def _parse_coordinates(self, data: Sequence[str]) -> List[Atom]:
        """Parse atomic coordinates from coordinate strings."""
        pattern = re.compile(
            r"^([a-z]{1,3})\s+"  # Element symbol (1-3 letters)
            r"(-?\d+\.?\d*)\s+"  # x coordinate
            r"(-?\d+\.?\d*)\s+"  # y coordinate
            r"(-?\d+\.?\d*)$",  # z coordinate
            re.IGNORECASE,
        )

        atoms = []
        for i, string in enumerate(data):
            if not (string := string.strip()):
                continue

            if not (match := pattern.match(string)):
                logger.error(f"Invalid coordinate format: {string}")
                raise ValueError(f"Invalid coordinate format: {string}")

            element, x, y, z = match.groups()
            if element not in self.default_radii:
                logger.error(f"Unknown element found: {element}")
                raise ValueError(f"Unknown element found: {element}")

            atom = Atom(
                element=element,
                position=Point3D(x=float(x), y=float(y), z=float(z)),
                index=i,
                radius=self.default_radii[element],
            )
            atoms.append(atom)

        if not atoms:
            logger.error("No valid coordinates found")
            raise ValueError("No valid coordinates found")

        return atoms

    def read_xyz(
        self, file_path: Union[str, Path], xyz_start: int = 2, validate: bool = False
    ) -> None:
        """Read molecular structure from XYZ file."""
        logger.debug(f"Reading XYZ file: {file_path}")

        file_path = Path(file_path)
        if not file_path.exists():
            logger.error(f"XYZ file not found: {file_path}")
            raise FileNotFoundError(f"XYZ file not found: {file_path}")

        try:
            lines = file_path.read_text().strip().split("\n")
            if len(lines) <= xyz_start:
                logger.error(f"XYZ file has no coordinate data after line {xyz_start}")
                raise ValueError(f"XYZ file has no coordinate data after line {xyz_start}")

            # Handle validation
            if validate and xyz_start > 0:
                expected_atoms = int(lines[0])

            # Store comment if available
            if xyz_start > 1:
                self.comment = lines[1].strip()

            # Parse coordinates and create atoms
            self.atoms = self._parse_coordinates(lines[xyz_start:])

            # Validate atom count if requested
            if validate and len(self.atoms) != expected_atoms:
                logger.error(
                    f"Number of coordinates doesn't match atom count in first line. "
                    f"Expected: {expected_atoms}, Found: {len(self.atoms)}"
                )
                raise ValueError(
                    f"Number of coordinates doesn't match atom count in first line. "
                    f"Expected: {expected_atoms}, Found: {len(self.atoms)}"
                )

            # Generate bonds
            self._generate_bonds()
            logger.info(f"Successfully read {len(self.atoms)} atoms from {file_path}")

        except Exception as e:
            logger.error(f"Error reading XYZ file: {e}", exc_info=True)
            raise

    def to_networkx(self) -> nx.Graph:
        """Convert to NetworkX graph."""
        logger.debug("Creating NetworkX graph")

        G = nx.Graph()

        # Add nodes with attributes
        for atom in self.atoms:
            G.add_node(atom.index, element=atom.element, xyz=(atom.x, atom.y, atom.z))

        # Add edges with attributes
        for bond in self.bonds:
            G.add_edge(bond.atom1.index, bond.atom2.index, length=bond.length)

        return G

    # Define dummy method to_plotly()
    def to_plotly(self, config: Optional[VisualizationConfig] = None) -> go.Figure:
        """Convert molecular structure to a Plotly figure.

        Args:
            config: Optional visualization configuration. If None, uses defaults.

        Returns:
            go.Figure: Plotly figure object with interactive molecular visualization
        """
        logger.debug("Creating Plotly figure")

        # Method 1: Using create_visualization helper function
        return create_visualization(self, config)

    def formula(self) -> str:
        """Return molecular formula in Hill notation."""
        if not self.atoms:
            return ""

        element_counts = Counter(atom.element for atom in self.atoms)
        formula_parts = []

        # Carbon and Hydrogen first
        if "C" in element_counts:
            count = element_counts.pop("C")
            formula_parts.append(f"C{count if count > 1 else ''}")

            if "H" in element_counts:
                count = element_counts.pop("H")
                formula_parts.append(f"H{count if count > 1 else ''}")

        # Add remaining elements alphabetically
        for element in sorted(element_counts):
            count = element_counts[element]
            formula_parts.append(f"{element}{count if count > 1 else ''}")

        return "".join(formula_parts)

    def __len__(self) -> int:
        """Return number of atoms."""
        return len(self.atoms)

    def __getitem__(self, index: int) -> Atom:
        """Get atom by index."""
        return self.atoms[index]

    def __iter__(self) -> Iterator[Tuple[str, Tuple[float, float, float]]]:
        """Create iterator over atoms returning (element, coordinates) tuples."""
        return ((atom.element, (atom.x, atom.y, atom.z)) for atom in self.atoms)
