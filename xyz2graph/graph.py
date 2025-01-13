"""Read, analyze, and visualize molecular structures from XYZ files.

Provides functionality for working with molecular graph structures through the MolGraph class,
including methods for manipulating element properties, reading XYZ files, and generating
visualizations.
"""

from collections import Counter
from dataclasses import dataclass, field
from itertools import compress
from pathlib import Path
from typing import (
    Dict,
    FrozenSet,
    List,
    Optional,
    Set,
    Union,
)

import networkx as nx
import numpy as np
import plotly.graph_objs as go
from numpy.typing import NDArray

from .constants import _DEFAULT_CPK_COLORS, _DEFAULT_RADII
from .logging import logger
from .molecule import Atom, Bond
from .visualization import VisualizationConfig, create_visualization


@dataclass
class MolGraph:
    """Represents a molecular graph structure.

    Provides functionality to read molecular structure data from XYZ files, analyze molecular
    geometry, and visualize the molecular structure using Plotly.

    Args:
        atoms (List[Atom], optional): List of Atom objects representing the molecular structure.
            Defaults to an empty list.
        bonds (List[Bond], optional): List of Bond objects representing molecular connectivity.
            Defaults to an empty list.
        comment (str, optional): Comment or description of the molecular structure.
            Defaults to an empty string.
        default_radii (Dict[str, float], optional): Dictionary mapping element symbols to their
            default atomic radii. Defaults to internal _DEFAULT_RADII.
        cpk_colors (Dict[str, str], optional): Dictionary mapping element symbols to their
            CPK color scheme colors. Defaults to internal _DEFAULT_CPK_COLORS.
        cpk_color_rest (str, optional): Default color for elements not specified in cpk_colors.
            Defaults to "pink".

    Attributes:
        elements (List[str]): List of element symbols in order of appearance.
        indices (List[int]): List of atom indices in order.
        x (List[float]): List of x coordinates for all atoms in order.
        y (List[float]): List of y coordinates for all atoms in order.
        z (List[float]): List of z coordinates for all atoms in order.
        xyz (NDArray[np.float64]): Nx3 array of atomic coordinates.
        atomic_radii (List[float]): List of atomic radii in atom order.
        bond_lengths (Dict[FrozenSet[int], float]): Mapping of atom index pairs to bond lengths.
        adj_list (Dict[int, Set[int]]): Dictionary mapping atom indices
            to sets of bonded atom indices.
        adj_matrix (NDArray[np.int_]): Square matrix representing molecular connectivity.
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
        """Get list of elements in the molecule.

        Returns:
            List[str]: List of element symbols in order of appearance.
        """
        return [atom.element for atom in self.atoms]

    @property
    def indices(self) -> List[int]:
        """Get list of atom indices.

        Returns:
            List[int]: List of atom indices in order.
        """
        return [atom.index for atom in self.atoms]

    @property
    def x(self) -> List[float]:
        """Get x coordinates of all atoms.

        Returns:
            List[float]: List of x coordinates in atom order.
        """
        return [atom.x for atom in self.atoms]

    @property
    def y(self) -> List[float]:
        """Get y coordinates of all atoms.

        Returns:
            List[float]: List of y coordinates in atom order.
        """
        return [atom.y for atom in self.atoms]

    @property
    def z(self) -> List[float]:
        """Get z coordinates of all atoms.

        Returns:
            List[float]: List of z coordinates in atom order.
        """
        return [atom.z for atom in self.atoms]

    @property
    def xyz(self) -> NDArray[np.float64]:
        """Get atomic coordinates as a numpy array.

        Returns:
            NDArray[np.float64]: Nx3 array of atomic coordinates.
        """
        return np.array([[atom.x, atom.y, atom.z] for atom in self.atoms])

    @property
    def atomic_radii(self) -> List[float]:
        """Get atomic radii for all atoms.

        Returns:
            List[float]: List of atomic radii in atom order.
        """
        return [atom.radius for atom in self.atoms]

    @property
    def bond_lengths(self) -> Dict[FrozenSet[int], float]:
        """Get dictionary of bond lengths.

        Returns:
            Dict[FrozenSet[int], float]: Mapping of atom index pairs to bond lengths.
        """
        return {bond.indices: bond.length for bond in self.bonds}

    @property
    def adj_list(self) -> Dict[int, Set[int]]:
        """Get adjacency list representation of molecular graph.

        Returns:
            Dict[int, Set[int]]: Dictionary mapping atom indices to sets of bonded atom indices.
        """
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

    def read_xyz(self, file_path: Union[str, Path]) -> None:
        """Read molecular structure from XYZ file.

        Args:
            file_path (Union[str, Path]): Path to XYZ format file.

        Raises:
            FileNotFoundError: If the specified file does not exist.
            ValueError: If file format is invalid or contains unknown elements.
        """
        file_path = Path(file_path)
        if not file_path.exists():
            raise FileNotFoundError(f"XYZ file not found: {file_path}")

        try:
            lines = file_path.read_text().splitlines()

            if len(lines) < 3:
                raise ValueError(
                    "XYZ file must contain at least 3 lines "
                    "(number of atoms, comment, and coordinates)"
                )
            try:
                n_atoms = int(lines[0])
            except (IndexError, ValueError) as err:
                raise ValueError("First line must be an integer (number of atoms)") from err

            self.comment = lines[1].strip()

            # Process coordinate lines
            coordinate_lines = list(filter(None, map(str.strip, lines[2:])))

            self.atoms = []
            for i, line in enumerate(coordinate_lines, start=0):
                parts = line.split()

                try:
                    element = parts[0]
                    if element not in self.default_radii:
                        raise ValueError(f"Unknown element symbol: {element}")

                    x, y, z = map(float, parts[1:4])

                    self.atoms.append(
                        Atom(
                            element=element,
                            x=x,
                            y=y,
                            z=z,
                            index=i,
                            radius=self.default_radii[element],
                        )
                    )
                except (IndexError, ValueError) as err:
                    raise ValueError(
                        f"Invalid format in line {i+3}, expected: element x y z"
                    ) from err

            if len(self.atoms) != n_atoms:
                logger.warning(
                    f"Number of atoms in file ({len(self.atoms)}) doesn't match the number "
                    f"specified in the first line ({n_atoms})"
                )

            self._generate_bonds()

        except Exception as e:
            logger.error(f"Error reading XYZ file: {e}")
            raise

    def _generate_bonds(self) -> None:
        """Generate bonds based on atomic distances.

        Updates the bonds list based on distance criteria between atoms.
        """
        self.bonds.clear()

        distances = self.distance_matrix()

        # Calculate bonding threshold matrix
        radii = np.array(self.atomic_radii)
        thresholds = (radii[:, np.newaxis] + radii) * 1.3

        # Find bonded atoms (creates adjacency matrix internally)
        adj_matrix = np.logical_and(0.1 < distances, thresholds > distances).astype(np.int_)

        # Create Bond objects from adjacency matrix
        for i, j in zip(*np.nonzero(adj_matrix)):
            if i < j:
                self.bonds.append(Bond(self.atoms[i], self.atoms[j]))

    def distance_matrix(self) -> NDArray[np.float64]:
        """Calculates the matrix of interatomic distances using optimized memory handling.

        Attempts to use fast vectorized calculation first, then falls back to a memory-efficient
        loop-based method if memory constraints are encountered.

        Returns:
            NDArray[np.float64]: A square matrix where entry (i,j) represents the
                Euclidean distance between atoms i and j in Angstroms.
        """
        try:
            distances = self.xyz[:, np.newaxis, :] - self.xyz
            return np.sqrt(np.einsum("ijk,ijk->ij", distances, distances))
        except MemoryError:
            # Fall back to loop-based method
            n_atoms = len(self.atoms)
            logger.info("Using memory-efficient method for distance calculation")
            distance_matrix = np.zeros((n_atoms, n_atoms), dtype=np.float64)
            for i in range(n_atoms):
                diff = self.xyz[i] - self.xyz
                distance_matrix[i] = np.sqrt(np.sum(diff * diff, axis=1))
            return distance_matrix

    def formula(self) -> str:
        """Generate molecular formula in Hill notation.

        Returns:
            str: Molecular formula in Hill notation.
        """
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

    def set_element_radius(self, element: str, radius: float) -> None:
        """Set the reference radius for a specific element.

        Args:
            element (str): Chemical element symbol.
            radius (float): New atomic radius value.
        """
        self.default_radii[element] = radius
        # Update radii for existing atoms of this element
        for atom in self.atoms:
            if atom.element == element:
                atom.radius = radius
        # Regenerate bonds with new radii
        self._generate_bonds()

    def remove(
        self,
        indices: Optional[List[int]] = None,
        elements: Optional[List[str]] = None,
        inplace: bool = False,
    ) -> Optional["MolGraph"]:
        """Remove atoms by indices and/or elements.

        Args:
            indices (List[int], optional): List of atom indices to remove.
            elements (List[str], optional): List of element symbols to remove.
            inplace (bool): If True, modify this instance. If False, return a new instance.

        Returns:
            Optional[MolGraph]: New MolGraph instance if inplace=False, None if inplace=True.

        Raises:
            IndexError: If any index is out of range.
            ValueError: If attempting to remove all atoms or if unknown elements specified.
        """
        mask = [True] * len(self.atoms)

        if indices is not None:
            if any(i < 0 or i >= len(self.atoms) for i in indices):
                raise IndexError("Atom index out of range")
            for idx in indices:
                mask[idx] = False

        if elements is not None:
            found_elements = {atom.element for atom in self.atoms if atom.element in elements}
            unused_elements = set(elements) - found_elements
            if unused_elements:
                logger.warning(
                    f"Element(s) not found: {', '.join(sorted(unused_elements))}. "
                    "Use proper case (e.g., 'H' not 'h')"
                )
            mask = [m and atom.element not in elements for m, atom in zip(mask, self.atoms)]

        if not any(mask):
            raise ValueError("Cannot remove all atoms from molecule")

        filtered_atoms = list(compress(self.atoms, mask))

        if inplace:
            self.atoms = filtered_atoms
            self._generate_bonds()
            return None

        new_mol = MolGraph()
        new_mol.atoms = filtered_atoms
        new_mol.comment = self.comment
        new_mol.default_radii = self.default_radii.copy()
        new_mol.cpk_colors = self.cpk_colors.copy()
        new_mol.cpk_color_rest = self.cpk_color_rest
        new_mol._generate_bonds()
        return new_mol

    def to_networkx(self) -> nx.Graph:
        """Convert molecular structure to NetworkX graph.

        Creates a NetworkX graph representation of the molecular structure where atoms are nodes
        and bonds are edges. Each node (atom) index corresponds to the atom's original index in
        the molecular structure and contains element type and 3D coordinates as attributes. Each
        edge (bond) connects two atoms using their indices and stores the bond length.

        Node attributes:
            - element (str): Chemical element symbol of the atom
            - xyz (tuple): 3D coordinates of the atom as (x, y, z) floats

        Edge attributes:
            - length (float): Bond length between the two connected atoms in Angstroms

        Returns:
            nx.Graph: NetworkX undirected graph where:
                - Nodes are atoms with their indices as node identifiers
                - Edges are bonds between atoms
                - Node and edge attributes contain chemical properties
        """
        logger.debug("Creating NetworkX graph")

        G = nx.Graph()

        # Add nodes with attributes
        for atom in self.atoms:
            G.add_node(atom.index, element=atom.element, xyz=(atom.x, atom.y, atom.z))

        # Add edges with attributes
        for bond in self.bonds:
            G.add_edge(bond.atom1.index, bond.atom2.index, length=bond.length)

        return G

    def to_plotly(self, config: Optional[VisualizationConfig] = None) -> go.Figure:
        """Convert molecular structure to Plotly figure.

        Args:
            config (VisualizationConfig, optional): Visualization configuration parameters.

        Returns:
            go.Figure: Interactive molecular visualization as Plotly figure.
        """
        logger.debug("Creating Plotly figure")

        return create_visualization(self, config)

    def __len__(self) -> int:
        """Get the number of atoms in the molecule.

        Returns:
            int: Total number of atoms.
        """
        return len(self.atoms)

    def __getitem__(self, index: int) -> Atom:
        """Get atom at specified index.

        Args:
            index (int): Zero-based index of the atom.

        Returns:
            Atom: Atom object at the specified index.

        Raises:
            IndexError: If index is out of range.
        """
        return self.atoms[index]

    def __repr__(self) -> str:
        """Creates a simplified string representation of the molecular graph.

        Generates a string representation showing the molecular formula in Hill notation
        along with basic structural information about the number of atoms and bonds.

        Returns:
            str: A string in the format "MolGraph(formula: num_atoms atoms, num_bonds bonds)"
                or "MolGraph(empty)" for an empty molecule.
        """
        if not self.atoms:
            return "MolGraph(empty)"

        return f"MolGraph({self.formula()}: {len(self.atoms)} atoms, {len(self.bonds)} bonds)"
