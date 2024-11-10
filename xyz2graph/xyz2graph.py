import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Counter, Dict, FrozenSet, Iterator, List, Sequence, Set, Tuple, Union

import networkx as nx
import numpy as np
import plotly.graph_objs as go
from numpy.typing import NDArray

from .logging import logger


DEFAULT_RADII: Dict[str, float] = {
    "Ac": 1.88,
    "Ag": 1.59,
    "Al": 1.35,
    "Am": 1.51,
    "As": 1.21,
    "Au": 1.50,
    "B": 0.83,
    "Ba": 1.34,
    "Be": 0.35,
    "Bi": 1.54,
    "Br": 1.21,
    "C": 0.68,
    "Ca": 0.99,
    "Cd": 1.69,
    "Ce": 1.83,
    "Cl": 0.99,
    "Co": 1.33,
    "Cr": 1.35,
    "Cs": 1.67,
    "Cu": 1.52,
    "D": 0.23,
    "Dy": 1.75,
    "Er": 1.73,
    "Eu": 1.99,
    "F": 0.64,
    "Fe": 1.34,
    "Ga": 1.22,
    "Gd": 1.79,
    "Ge": 1.17,
    "H": 0.23,
    "Hf": 1.57,
    "Hg": 1.70,
    "Ho": 1.74,
    "I": 1.40,
    "In": 1.63,
    "Ir": 1.32,
    "K": 1.33,
    "La": 1.87,
    "Li": 0.68,
    "Lu": 1.72,
    "Mg": 1.10,
    "Mn": 1.35,
    "Mo": 1.47,
    "N": 0.68,
    "Na": 0.97,
    "Nb": 1.48,
    "Nd": 1.81,
    "Ni": 1.50,
    "Np": 1.55,
    "O": 0.68,
    "Os": 1.37,
    "P": 1.05,
    "Pa": 1.61,
    "Pb": 1.54,
    "Pd": 1.50,
    "Pm": 1.80,
    "Po": 1.68,
    "Pr": 1.82,
    "Pt": 1.50,
    "Pu": 1.53,
    "Ra": 1.90,
    "Rb": 1.47,
    "Re": 1.35,
    "Rh": 1.45,
    "Ru": 1.40,
    "S": 1.02,
    "Sb": 1.46,
    "Sc": 1.44,
    "Se": 1.22,
    "Si": 1.20,
    "Sm": 1.80,
    "Sn": 1.46,
    "Sr": 1.12,
    "Ta": 1.43,
    "Tb": 1.76,
    "Tc": 1.35,
    "Te": 1.47,
    "Th": 1.79,
    "Ti": 1.47,
    "Tl": 1.55,
    "Tm": 1.72,
    "U": 1.58,
    "V": 1.33,
    "W": 1.37,
    "Y": 1.78,
    "Yb": 1.94,
    "Zn": 1.45,
    "Zr": 1.56,
}

DEFAULT_CPK_COLORS: Dict[str, str] = {
    "Ar": "cyan",
    "B": "salmon",
    "Ba": "darkgreen",
    "Be": "darkgreen",
    "Br": "darkred",
    "C": "black",
    "Ca": "darkgreen",
    "Cl": "green",
    "Cs": "violet",
    "F": "green",
    "Fe": "darkorange",
    "Fr": "violet",
    "H": "white",
    "He": "cyan",
    "I": "darkviolet",
    "K": "violet",
    "Kr": "cyan",
    "Li": "violet",
    "Mg": "darkgreen",
    "N": "blue",
    "Na": "violet",
    "Ne": "cyan",
    "O": "red",
    "P": "orange",
    "Ra": "darkgreen",
    "Rb": "violet",
    "S": "yellow",
    "Sr": "darkgreen",
    "Ti": "gray",
    "Xe": "cyan",
}


@dataclass
class MolGraph:
    """
    Represents a molecular graph structure from XYZ file format.

    This class provides functionality to read molecular structure data from XYZ files,
    analyze molecular geometry, and visualize the molecular structure using Plotly.
    Colors and atomic radii can be customized per instance.

    Attributes:
        elements: List of atomic elements in the molecule
        x: List of x-coordinates for each atom
        y: List of y-coordinates for each atom
        z: List of z-coordinates for each atom
        adj_list: Dictionary mapping atom indices to their connected neighbors
        atomic_radii: List of atomic radii for each atom
        bond_lengths: Dictionary mapping pairs of connected atoms to their bond lengths
        adj_matrix: Numpy array representing the adjacency matrix of the molecular graph
        default_radii: Dictionary of reference atomic radii for each element
        cpk_colors: Dictionary of CPK colors for each element
        cpk_color_rest: Default color for elements not in cpk_colors
    """

    elements: List[str] = field(default_factory=list)
    x: List[float] = field(default_factory=list)
    y: List[float] = field(default_factory=list)
    z: List[float] = field(default_factory=list)
    adj_list: Dict[int, Set[int]] = field(default_factory=dict)
    atomic_radii: List[float] = field(default_factory=list)
    bond_lengths: Dict[FrozenSet[int], float] = field(default_factory=dict)
    adj_matrix: Union[NDArray[np.int_], None] = field(default=None)

    # Customizable parameters with defaults
    default_radii: Dict[str, float] = field(default_factory=lambda: dict(DEFAULT_RADII))
    cpk_colors: Dict[str, str] = field(default_factory=lambda: dict(DEFAULT_CPK_COLORS))
    cpk_color_rest: str = field(default="pink")

    def set_element_radius(self, element: str, radius: float) -> None:
        """
        Set the reference radius for a specific element and update the adjacency list.

        Args:
            element: Chemical element symbol
            radius: New radius value in Angstroms
        """
        self.default_radii[element] = radius
        if element in self.elements:
            self.atomic_radii = [
                self.default_radii[element] for element in self.elements
            ]
            self._generate_adjacency_list()

    def set_element_color(self, element: str, color: str) -> None:
        """
        Set the CPK color for a specific element.

        Args:
            element: Chemical element symbol
            color: Color name or code (any format accepted by Plotly)
        """
        self.cpk_colors[element] = color

    def set_default_color(self, color: str) -> None:
        """
        Set the default color for elements not in cpk_colors.

        Args:
            color: Color name or code (any format accepted by Plotly)
        """
        self.cpk_color_rest = color

    def _parse_coordinates(
        self, data: Sequence[str]
    ) -> Tuple[List[str], List[float], List[float], List[float]]:
        """
        Parse atomic coordinates from a list of coordinate strings.

        Args:
            data: List of strings, each containing element and coordinates in format:
                        'element x y z' (e.g., 'H 0.0 0.0 0.0')

        Returns:
            Tuple containing (elements, x_coords, y_coords, z_coords)

        Raises:
            ValueError: If coordinate format is invalid or contains unknown elements
        """
        elements, xs, ys, zs = [], [], [], []

        pattern = re.compile(
            r"^([a-z]{1,3})\s+"  # Element symbol (1-3 letters)
            r"(-?\d+\.?\d*)\s+"  # x coordinate (integer or float)
            r"(-?\d+\.?\d*)\s+"  # y coordinate (integer or float)
            r"(-?\d+\.?\d*)$",  # z coordinate (integer or float)
            re.IGNORECASE,
        )

        for string in data:
            if not (string := string.strip()):
                continue

            if not (match := pattern.match(string)):
                logger.error(f"Invalid coordinate format: {string}")
                raise ValueError(f"Invalid coordinate format: {string}")

            element, x, y, z = match.groups()
            if element not in self.default_radii:
                logger.error(f"Unknown element found: {element}")
                raise ValueError(f"Unknown element found: {element}")

            elements.append(element)
            xs.append(float(x))
            ys.append(float(y))
            zs.append(float(z))

        if not elements:
            logger.error("No valid coordinates found")
            raise ValueError("No valid coordinates found")

        return elements, xs, ys, zs

    def read_xyz(
        self, file_path: Union[str, Path], xyz_start: int = 2, validate: bool = False
    ) -> None:
        """
        Reads molecular structure data from an XYZ file.

        The XYZ file format specification:
        - Line 0: Number of atoms (integer)
        - Line 1: Comment or title (can be empty)
        - Line 2 onwards: Atomic coordinates in format:
            element  x  y  z
        where:
        - element: Chemical element symbol (e.g., H, He, Li)
        - x, y, z: Cartesian coordinates in Angstroms (float)

        Example of XYZ file:
        3
        Water molecule
        O  0.000000  0.000000  0.000000
        H  0.758602  0.000000  0.504284
        H  0.758602  0.000000 -0.504284

        Args:
            file_path: Path to the XYZ file
            xyz_start: Line number where XYZ coordinate data starts (default=2)
            validate: If True, validates that the number of coordinates matches
                    the atom count in the first line. Note: validation is automatically
                    disabled if xyz_start=0

        Raises:
            FileNotFoundError: If the specified file does not exist
            ValueError: If the file format is invalid or contains unknown elements
        """
        logger.debug(f"Reading XYZ file: {file_path}")

        # Handle validation conflict
        if validate and xyz_start < 1:
            logger.warning(
                "Atom count validation disabled:"
                " cannot validate when reading from first line "
                "(xyz_start=0)"
            )
            validate = False

        file_path = Path(file_path)
        if not file_path.exists():
            logger.error(f"XYZ file not found: {file_path}")
            raise FileNotFoundError(f"XYZ file not found: {file_path}")

        try:
            lines = file_path.read_text().strip().split("\n")
            if len(lines) <= xyz_start:
                logger.error(f"XYZ file has no coordinate data after line {xyz_start}")
                raise ValueError(
                    f"XYZ file has no coordinate data after line {xyz_start}"
                )

            # Get expected atom count if validation is requested
            expected_atoms = None
            if validate:
                try:
                    expected_atoms = int(lines[0])
                    logger.debug(f"Expected number of atoms: {expected_atoms}")
                except (IndexError, ValueError) as err:
                    logger.error("Could not read atom count from first line")
                    raise ValueError(
                        "Could not read atom count from first line"
                    ) from err

            coordinates = self._parse_coordinates(lines[xyz_start:])
            self.elements, self.x, self.y, self.z = coordinates

            # Validate number of atoms if requested
            if validate and len(self.elements) != expected_atoms:
                logger.error(
                    f"Number of coordinates doesn't match atom count in first line. "
                    f"Expected: {expected_atoms}, Found: {len(self.elements)}"
                )
                raise ValueError(
                    f"Number of coordinates doesn't match atom count in first line. "
                    f"Expected: {expected_atoms}, Found: {len(self.elements)}"
                )

            logger.info(
                f"Successfully read {len(self.elements)} atoms from {file_path}"
            )

            self.atomic_radii = [
                self.default_radii[element] for element in self.elements
            ]
            self._generate_adjacency_list()

        except Exception as e:
            logger.error(f"Error reading XYZ file: {e}", exc_info=True)
            raise

    def to_plotly(self) -> go.Figure:
        """
        Creates a Plotly figure for 3D visualization of the molecule.

        Returns:
            go.Figure: Interactive 3D visualization of the molecular structure
        """

        logger.debug("Creating Plotly figure")

        def atom_trace() -> go.Scatter3d:
            """Creates an atom trace for the plot."""
            colors = [
                self.cpk_colors.get(element, self.cpk_color_rest)
                for element in self.elements
            ]
            markers = dict(
                color=colors,
                line=dict(color="lightgray", width=2),
                size=7,
                symbol="circle",
                opacity=0.8,
            )
            return go.Scatter3d(
                x=self.x,
                y=self.y,
                z=self.z,
                mode="markers",
                marker=markers,
                text=self.elements,
            )

        def bond_trace() -> go.Scatter3d:
            """Creates a bond trace for the plot."""
            trace = go.Scatter3d(
                x=[],
                y=[],
                z=[],
                hoverinfo="none",
                mode="lines",
                marker=dict(color="grey", size=7, opacity=1),
            )
            for i, j in self.edges():
                trace["x"] += (self.x[i], self.x[j], None)
                trace["y"] += (self.y[i], self.y[j], None)
                trace["z"] += (self.z[i], self.z[j], None)
            return trace

        # Create annotations
        def create_annotations(
            elements: bool = True, indices: bool = False, bond_lengths: bool = False
        ) -> List[dict]:
            annotations = []

            if elements or indices:
                for idx, (element, (x, y, z)) in enumerate(self):
                    text = element if elements else str(idx)
                    annotations.append(
                        dict(text=text, x=x, y=y, z=z, showarrow=False, yshift=15)
                    )

            if bond_lengths:
                for (i, j), length in self.bond_lengths.items():
                    x = (self.x[i] + self.x[j]) / 2
                    y = (self.y[i] + self.y[j]) / 2
                    z = (self.z[i] + self.z[j]) / 2
                    annotations.append(
                        dict(
                            text=round(length, 2),
                            x=x,
                            y=y,
                            z=z,
                            showarrow=False,
                            yshift=15,
                            font=dict(color="steelblue"),
                        )
                    )

            return annotations

        # Create update menu
        buttons = [
            dict(
                label=" Elements",
                method="relayout",
                args=[{"scene.annotations": create_annotations(elements=True)}],
            ),
            dict(
                label=" Elements & Bond Lengths",
                method="relayout",
                args=[
                    {
                        "scene.annotations": create_annotations(
                            elements=True, bond_lengths=True
                        )
                    }
                ],
            ),
            dict(
                label="Indices",
                method="relayout",
                args=[{"scene.annotations": create_annotations(indices=True)}],
            ),
            dict(
                label="Indices & Bond Lengths",
                method="relayout",
                args=[
                    {
                        "scene.annotations": create_annotations(
                            indices=True, bond_lengths=True
                        )
                    }
                ],
            ),
            dict(
                label="Bond Lengths",
                method="relayout",
                args=[{"scene.annotations": create_annotations(bond_lengths=True)}],
            ),
            dict(
                label="Hide All",
                method="relayout",
                args=[{"scene.annotations": []}],
            ),
        ]

        updatemenus = [
            dict(
                buttons=buttons,
                direction="down",
                xanchor="left",
                yanchor="top",
            )
        ]

        # Create figure
        data = [atom_trace(), bond_trace()]
        axis_params = dict(
            showgrid=False,
            showbackground=False,
            showticklabels=False,
            zeroline=False,
            titlefont=dict(color="white"),
        )
        layout = dict(
            scene=dict(
                xaxis=axis_params,
                yaxis=axis_params,
                zaxis=axis_params,
                annotations=create_annotations(),
            ),
            margin=dict(r=0, l=0, b=0, t=0),
            showlegend=False,
            updatemenus=updatemenus,
        )
        logger.debug("Plotly figure created successfully")
        return go.Figure(data=data, layout=layout)

    def to_networkx(self) -> nx.Graph:
        """
        Converts the molecular structure to a NetworkX graph.

        Returns:
            nx.Graph: NetworkX graph representation with node and edge attributes
        """
        logger.debug("Creating NetworkX graph")

        G = nx.Graph(self.adj_list)

        # Add node attributes
        node_attrs = {
            num: {"element": element, "xyz": (self.x[num], self.y[num], self.z[num])}
            for num, element in enumerate(self.elements)
        }
        nx.set_node_attributes(G, node_attrs)

        # Add edge attributes
        edge_attrs = {
            edge: {"length": length} for edge, length in self.bond_lengths.items()
        }
        nx.set_edge_attributes(G, edge_attrs)
        logger.debug("NetworkX graph created successfully")

        return G

    def _generate_adjacency_list(self) -> None:
        """
        Generates the adjacency list and matrix based on atomic positions and radii.

        This method uses atomic positions and reference radii to determine bonding
        between atoms. Two atoms are considered bonded if their distance is less
        than 1.3 times the sum of their atomic radii.
        """
        logger.debug("Generating adjacency list")

        try:
            xyz = np.stack((self.x, self.y, self.z), axis=-1)
            distances = xyz[:, np.newaxis, :] - xyz
            distances = np.sqrt(np.einsum("ijk,ijk->ij", distances, distances))

            atomic_radii = np.array(self.atomic_radii)
            distance_bond = (atomic_radii[:, np.newaxis] + atomic_radii) * 1.3

            self.adj_matrix = np.logical_and(
                0.1 < distances, distance_bond > distances
            ).astype(int)

            # Clear existing adjacency data
            self.adj_list.clear()
            self.bond_lengths.clear()

            # Generate new adjacency data
            for i, j in zip(*np.nonzero(self.adj_matrix)):
                self.adj_list.setdefault(i, set()).add(j)
                self.adj_list.setdefault(j, set()).add(i)
                self.bond_lengths[frozenset([i, j])] = round(distances[i, j], 5)
            logger.debug("Adjacency list generated successfully")

        except Exception as e:
            logger.error(f"Error generating adjacency list: {e}", exc_info=True)
            raise

    def edges(self) -> Iterator[Tuple[int, int]]:
        """
        Creates an iterator with all graph edges.

        Returns:
            Iterator[Tuple[int, int]]: Iterator yielding pairs of connected atom indices
        """
        edges = set()
        for node, neighbours in self.adj_list.items():
            for neighbour in neighbours:
                edge = frozenset([node, neighbour])
                if edge in edges:
                    continue
                edges.add(edge)
                yield node, neighbour

    def formula(self) -> str:
        """
        Returns the molecular formula in standard Hill notation.

        The Hill system orders elements as:
        1. Carbon first (if present)
        2. Hydrogen second (if carbon is present)
        3. All other elements alphabetically

        Returns:
            str: Molecular formula string
        """
        if not self.elements:
            return ""

        element_counts = Counter(self.elements)

        # Build formula in Hill notation
        formula_parts = []

        # Carbon and Hydrogen first if carbon is present
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
        """
        Returns the number of atoms in the molecule.

        Returns:
            int: Number of atoms
        """
        return len(self.elements)

    def __getitem__(self, position: int) -> Tuple[str, Tuple[float, float, float]]:
        """
        Returns the element and coordinates for the atom at the given position.

        Args:
            position: Index of the atom

        Returns:
            Tuple[str, Tuple[float, float, float]]: Tuple containing the element
            symbol and its (x, y, z) coordinates

        Raises:
            IndexError: If position is out of range
        """
        if not 0 <= position < len(self):
            raise IndexError("Atom index out of range")

        return self.elements[position], (
            self.x[position],
            self.y[position],
            self.z[position],
        )

    def __iter__(self) -> Iterator[Tuple[str, Tuple[float, float, float]]]:
        """
        Creates an iterator over all atoms in the molecule.

        Returns:
            Iterator[Tuple[str, Tuple[float, float, float]]]: Iterator yielding
            tuples of (element, coordinates) for each atom
        """
        return (self[i] for i in range(len(self)))

    def __repr__(self) -> str:
        """
        Creates a simplified string representation of the molecular graph.

        Returns:
            str: A string showing molecular formula and basic structural information
        """
        if not self.elements:
            return "MolGraph(empty)"

        element_counts = Counter(self.elements)

        # Build formula prioritizing C and H
        formula_parts = []
        for element in ["C", "H"]:
            if count := element_counts.pop(element, 0):
                formula_parts.append(f"{element}{count if count > 1 else ''}")

        # Add remaining elements alphabetically
        for element in sorted(element_counts):
            count = element_counts[element]
            formula_parts.append(f"{element}{count if count > 1 else ''}")

        formula = "".join(formula_parts)
        num_bonds = sum(len(neighbors) for neighbors in self.adj_list.values()) // 2

        return f"MolGraph({formula}: {len(self.elements)} atoms, {num_bonds} bonds)"
