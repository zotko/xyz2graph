"""Read, analyze, and visualize molecular structures from XYZ files.

It includes a `MolGraph` class that represents a molecular graph structure, with methods for
setting element radii and colors, filtering atoms, reading XYZ files, generating adjacency lists,
and converting to Plotly and NetworkX representations.

Classes:
    PlotConfig: TypedDict for configuring molecular plot visualization.
    MolGraph: Represents a molecular graph structure with methods for reading, analyzing, and
        visualizing molecular data.

Constants:
    DEFAULT_RADII: Default atomic radii for various elements.
    DEFAULT_CPK_COLORS: Default CPK colors for various elements.
    DEFAULT_PLOT_CONFIG: Default configuration for molecular plot visualization.

Functions:
    set_element_radius: Sets the reference radius for a specific element and updates the adjacency
        list.
    set_element_color: Sets the CPK color for a specific element.
    set_default_color: Sets the default color for elements not in cpk_colors.
    filter: Filters out atoms by indices and/or elements.
    _parse_coordinates: Parses atomic coordinates from a list of coordinate strings.
    read_xyz: Reads molecular structure data from an XYZ file.
    _create_atom_trace: Creates a Plotly trace for atoms in the molecule.
    _create_bond_trace: Creates a Plotly trace for bonds in the molecule.
    _create_annotations: Creates annotations for atom and bond labels in the plot.
    _create_menu_buttons: Creates menu buttons for toggling different label displays.
    to_plotly: Creates a Plotly figure for 3D visualization of the molecule.
    to_networkx: Converts the molecular structure to a NetworkX graph.
    _generate_adjacency_list: Generates the adjacency list and matrix based on atomic positions
        and radii.
    bonds: Creates an iterator with all graph bonds.
    edges: Creates an iterator with all graph bonds (alias for bonds method).
    formula: Returns the molecular formula in standard Hill notation.
    __len__: Returns the number of atoms in the molecule.
    __getitem__: Returns the element and coordinates for the atom at the given position.
    __iter__: Creates an iterator over all atoms in the molecule.
    __repr__: Creates a simplified string representation of the molecular graph.
"""

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import (
    Counter,
    Dict,
    FrozenSet,
    Iterator,
    List,
    Optional,
    Sequence,
    Set,
    Tuple,
    TypedDict,
    Union,
)

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


class PlotConfig(TypedDict, total=False):
    """Configuration dictionary for molecular plot visualization.

    Attributes:
        atom_size: Size of atom markers in the plot. Default is 7.
        atom_opacity: Opacity of atom markers (0.0 to 1.0). Default is 0.8.
        atom_border_color: Color string for atom marker borders. Default is "lightgray".
        atom_border_width: Width of atom marker borders in pixels. Default is 2.
        bond_color: Color string for bonds between atoms. Default is "grey".
        bond_width: Width of bond lines in pixels. Default is 2.
        show_grid: Whether to display the plot grid. Default is False.
        label_offset: Vertical offset for labels in pixels. Default is 15.
        bond_label_color: Color string for bond length labels. Default is "steelblue".
    """

    atom_size: int
    atom_opacity: float
    atom_border_color: str
    atom_border_width: int
    bond_color: str
    bond_width: int
    show_grid: bool
    label_offset: int
    bond_label_color: str


DEFAULT_PLOT_CONFIG: PlotConfig = {
    "atom_size": 7,
    "atom_opacity": 0.8,
    "atom_border_color": "lightgray",
    "atom_border_width": 2,
    "bond_color": "grey",
    "bond_width": 2,
    "show_grid": False,
    "label_offset": 15,
    "bond_label_color": "steelblue",
}


@dataclass
class MolGraph:
    """Represents a molecular graph structure.

    This class provides functionality to read molecular structure data from XYZ files,
    analyze molecular geometry, and visualize the molecular structure using Plotly.

    Attributes:
        elements (List[str]): List of atomic elements in the molecule.
        x (List[float]): List of x-coordinates for each atom.
        y (List[float]): List of y-coordinates for each atom.
        z (List[float]): List of z-coordinates for each atom.
        comment (str): Optional comment or title from the XYZ file.
        adj_list (Dict[int, Set[int]]): Dictionary mapping atom indices to their
            connected neighbors.
        atomic_radii (List[float]): List of atomic radii for each atom.
        bond_lengths (Dict[FrozenSet[int], float]): Dictionary mapping pairs of
            connected atoms to their bond lengths.
        adj_matrix (NDArray[np.int_] | None): Numpy array representing the adjacency
            matrix of the molecular graph.
        indices (List[int]): List of original indices for each atom, preserved during filtering
            operations. Maps current atoms back to their positions in the original structure.

    Examples:
        >>> mg = MolGraph()
        >>> mg.read_xyz('molecule.xyz')
        >>> print(mg.formula())
        'C2H6O'
    """

    elements: List[str] = field(default_factory=list)
    x: List[float] = field(default_factory=list)
    y: List[float] = field(default_factory=list)
    z: List[float] = field(default_factory=list)
    comment: str = field(default="")
    adj_list: Dict[int, Set[int]] = field(default_factory=dict)
    atomic_radii: List[float] = field(default_factory=list)
    bond_lengths: Dict[FrozenSet[int], float] = field(default_factory=dict)
    adj_matrix: Union[NDArray[np.int_], None] = field(default=None)
    indices: List[int] = field(default_factory=list)

    # Customizable parameters with defaults
    default_radii: Dict[str, float] = field(default_factory=lambda: dict(DEFAULT_RADII))
    cpk_colors: Dict[str, str] = field(default_factory=lambda: dict(DEFAULT_CPK_COLORS))
    cpk_color_rest: str = field(default="pink")

    def set_element_radius(self, element: str, radius: float) -> None:
        """Set the reference radius for a specific element and update the adjacency list.

        Updates the reference radius for the specified chemical element
        and recalculates the molecular adjacency list based on the new radius value.

        Args:
            element (str): Chemical element symbol
            radius (float): New radius value in Angstroms

        Examples:
            >>> mg.set_element_radius('C', 0.75)
        """
        self.default_radii[element] = radius
        if element in self.elements:
            self.atomic_radii = [self.default_radii[element] for element in self.elements]
            self._generate_adjacency_list()

    def set_element_color(self, element: str, color: str) -> None:
        """Set the CPK color for a specific element.

        Assigns a custom CPK (Corey-Pauling-Koltun) color to the specified chemical element
        for use in molecular visualizations.

        Args:
            element (str): Chemical element symbol
            color (str): Color name or code (any format accepted by Plotly)

        Examples:
            >>> mg.set_element_color('O', 'green')
            >>> mg.set_element_color('C', '#00FF00')
        """
        self.cpk_colors[element] = color

    def set_default_color(self, color: str) -> None:
        """Set the default color for elements not in the standard CPK color scheme.

        Specifies the fallback color to use for any chemical elements
        that don't have a standard CPK color assignment.

        Args:
            color (str): Color name or code (any format accepted by Plotly)

        Examples:
            >>> mg.set_default_color('gray')
            >>> mg.set_default_color('#A0A0A0')
        """
        self.cpk_color_rest = color

    def filter(
        self,
        indices: Optional[List[int]] = None,
        elements: Optional[List[str]] = None,
        inplace: bool = False,
    ) -> Optional["MolGraph"]:
        """Filter out atoms by indices and/or elements.

        Filter atoms from the molecular structure based on their indices or element types.
        Can either modify the existing structure or create a new filtered copy.

        Args:
            indices (List[int], optional): List of atom indices to filter out.
                If None, no indices are removed.
            elements (List[str], optional): List of element symbols to filter out.
                If None, no elements are removed.
            inplace (bool, default=False): If True, modifies the current object.
                If False, returns a new object.

        Returns:
            Optional[MolGraph]: A new filtered MolGraph object if inplace=False,
                None if inplace=True.

        Raises:
            ValueError: If filtering would remove all atoms.
            IndexError: If any index is out of range.

        Examples:
            >>> mg_without_h = mg.filter(elements=['H']) # Filter out all hydrogen atoms
            >>> mg_trimmed = mg.filter(indices=[0, 1, 2]) # Filter out specific atoms by index
            >>> mg_modified = mg.filter(indices=[0, 1, 2], elements=['H'])
            >>> mg.filter(elements=['H'], inplace=True) # Modify in place
        """
        from itertools import compress, starmap
        from operator import and_

        # Create mask for all atoms (True means keep)
        mask = [True] * len(self)

        # Apply indices filter if specified
        if indices is not None:
            # Validate indices
            if any(i < 0 or i >= len(self) for i in indices):
                raise IndexError("Atom index out of range")
            for idx in indices:
                mask[idx] = False

        # Apply elements filter if specified
        if elements is not None:
            # Track which filter elements are actually used
            found_elements = set()
            for elem in self.elements:
                if elem in elements:
                    found_elements.add(elem)

            # Warn about unused filter elements
            unused_elements = set(elements) - found_elements
            if unused_elements:
                logger.warning(
                    f"Element(s) not found: {', '.join(sorted(unused_elements))}. "
                    "Use proper case (e.g., 'H' not 'h')"
                )

        # Update mask to remove specified elements
        element_mask = [elem not in (elements or []) for elem in self.elements]
        mask = list(starmap(and_, zip(mask, element_mask)))

        # Check if any atoms remain
        if not any(mask):
            raise ValueError("Cannot filter out all atoms from molecule")

        # If modifying in place
        if inplace:
            # Store filtered data using compress
            self.elements = list(compress(self.elements, mask))
            self.x = list(compress(self.x, mask))
            self.y = list(compress(self.y, mask))
            self.z = list(compress(self.z, mask))
            self.atomic_radii = list(compress(self.atomic_radii, mask))
            self.indices = list(compress(self.indices, mask))

            # Regenerate adjacency information
            self._generate_adjacency_list()
            return None

        # Create new instance
        new_mg = MolGraph()
        new_mg.elements = list(compress(self.elements, mask))
        new_mg.x = list(compress(self.x, mask))
        new_mg.y = list(compress(self.y, mask))
        new_mg.z = list(compress(self.z, mask))
        new_mg.atomic_radii = list(compress(self.atomic_radii, mask))
        new_mg.indices = list(compress(self.indices, mask))
        new_mg.comment = self.comment

        # Copy customization settings
        new_mg.default_radii = self.default_radii.copy()
        new_mg.cpk_colors = self.cpk_colors.copy()
        new_mg.cpk_color_rest = self.cpk_color_rest

        # Generate adjacency information for new molecule
        new_mg._generate_adjacency_list()
        return new_mg

    def _parse_coordinates(
        self, data: Sequence[str]
    ) -> Tuple[List[str], List[float], List[float], List[float]]:
        """Parse atomic coordinates from a coordinate string list.

        Takes a list of coordinate strings in 'element x y z' format and parses them into
        separate lists of elements and their corresponding x, y, z coordinates. Each coordinate
        string should contain an element symbol followed by three space-separated numeric values
        representing the atomic position in 3D space.

        Args:
            data (Sequence[str]): List of coordinate strings, where each string contains an
            element symbol and coordinates in the format 'element x y z' (e.g., 'H 0.0 0.0 0.0').

        Returns:
            Tuple[List[str], List[float], List[float], List[float]]: A tuple containing four lists:
            elements, x_coordinates, y_coordinates, and z_coordinates.

        Raises:
            ValueError: If any coordinate string has invalid format or contains unknown element
            symbols not present in the default radii dictionary.

        Examples:
            >>> mg = MolGraph()
            >>> coords = ['H 0.0 0.0 0.0', 'O 0.0 0.0 1.0']
            >>> elements, x, y, z = mg._parse_coordinates(coords)
            >>> print(elements)
            ['H', 'O']
            >>> print(x)
            [0.0, 0.0]
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
        """Reads molecular structure data from an XYZ file.

        The file should have a header with the number of atoms and optional comment, followed by
        atomic coordinates in Angstroms with element symbols.

        Args:
            file_path (Union[str, Path]): Path to the XYZ file
            xyz_start (int, optional): Line number where XYZ coordinate data starts. Defaults to 2.
            validate (bool, optional, default=False): If True, validates that the number
                of coordinates matches the atom count in the first line. Note: validation is
                automatically disabled if xyz_start=0. Defaults to False.

        Raises:
            FileNotFoundError: If the specified file does not exist
            ValueError: If the file format is invalid or contains unknown elements

        Notes:
            The XYZ file format is a standard chemical file format that contains:

                - First line: Number of atoms
                - Second line: Comment or molecule name (optional)
                - Remaining lines: Atomic coordinates in format "element x y z"

            The function supports both standard XYZ files (xyz_start=2) and raw coordinate data
            (xyz_start=0). When validation is enabled, it checks that the number of coordinates
            matches the atom count specified in the first line.

        Examples:
            >>> mg = MolGraph()
            >>> mg.read_xyz("molecule.xyz")  # Read standard XYZ file (starts at line 2)
            >>> mg.read_xyz("coords.xyz", xyz_start=0)  # Read XYZ data starting from first line
            >>> mg.read_xyz("molecule.xyz", validate=True)  # Enable validation of atom count
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
                raise ValueError(f"XYZ file has no coordinate data after line {xyz_start}")

            # Get expected atom count if validation is requested
            expected_atoms = None
            if validate:
                try:
                    expected_atoms = int(lines[0])
                    logger.debug(f"Expected number of atoms: {expected_atoms}")
                except (IndexError, ValueError) as err:
                    logger.error("Could not read atom count from first line")
                    raise ValueError("Could not read atom count from first line") from err

            # Store comment if available (line 1 in standard XYZ format)
            if xyz_start > 1:
                self.comment = lines[1].strip()
                logger.debug(f"Read comment: {self.comment}")

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

            logger.info(f"Successfully read {len(self.elements)} atoms from {file_path}")

            self.atomic_radii = [self.default_radii[element] for element in self.elements]
            self._generate_adjacency_list()

            # Store only essential geometric and reference data
            self._original_state = {
                "elements": self.elements.copy(),
                "x": self.x.copy(),
                "y": self.y.copy(),
                "z": self.z.copy(),
                "indices": self.indices.copy(),
            }

        except Exception as e:
            logger.error(f"Error reading XYZ file: {e}", exc_info=True)
            raise

    def _create_atom_trace(self, plot_config: PlotConfig) -> go.Scatter3d:
        """Creates a Plotly trace for atoms in the molecule."""
        colors = [self.cpk_colors.get(element, self.cpk_color_rest) for element in self.elements]

        hover_text = [
            f"{i}. {elem} ({x:.2f}, {y:.2f}, {z:.2f})"
            for i, (elem, x, y, z) in enumerate(zip(self.elements, self.x, self.y, self.z))
        ]

        return go.Scatter3d(
            x=self.x,
            y=self.y,
            z=self.z,
            mode="markers",
            marker=dict(
                color=colors,
                line=dict(
                    color=plot_config["atom_border_color"],
                    width=plot_config["atom_border_width"],
                ),
                size=plot_config["atom_size"],
                opacity=plot_config["atom_opacity"],
                symbol="circle",
            ),
            text=self.elements,
            hovertext=hover_text,
            hoverinfo="text",
            name="atoms",
        )

    def _create_bond_trace(self, plot_config: PlotConfig) -> go.Scatter3d:
        """Creates a Plotly trace for bonds in the molecule."""
        # Create line segments for each bond (use None to separate segments)
        xs, ys, zs = [], [], []
        for i, j in self.edges():
            xs.extend([self.x[i], self.x[j], None])
            ys.extend([self.y[i], self.y[j], None])
            zs.extend([self.z[i], self.z[j], None])

        return go.Scatter3d(
            x=xs,
            y=ys,
            z=zs,
            mode="lines",
            line=dict(color=plot_config["bond_color"], width=plot_config["bond_width"]),
            hoverinfo="none",
            name="bonds",
        )

    def _create_annotations(
        self,
        plot_config: PlotConfig,
        atom_label_type: Optional[str] = None,
        show_distances: bool = False,
    ) -> List[dict]:
        """Create annotations for atom and bond labels in the molecular plot.

        Creates Plotly annotation dictionaries for displaying atom labels (either element
        symbols or indices) and/or bond length labels. The annotations are positioned above
        each atom or bond with configurable vertical offset.

        Args:
            plot_config (PlotConfig): Plot configuration dictionary containing visual
                parameters like label offsets and colors.
            atom_label_type (Optional[str]): Type of atom labels to show. Must be one of:
                None, "element", or "index". None means no atom labels.
            show_distances (bool): Whether to show bond length labels centered on each bond.

        Returns:
            List[dict]: List of Plotly annotation dictionaries containing text, position,
                and style information for atom and bond labels.

        Raises:
            ValueError: If atom_label_type is not one of: None, "element", or "index".
        """
        labels = []

        if atom_label_type is not None:
            if atom_label_type not in ("element", "index"):
                raise ValueError('atom_label_type must be None, "element", or "index"')

            for idx, (element, (x, y, z)) in enumerate(self):
                text = element if atom_label_type == "element" else str(idx)
                labels.append(
                    dict(
                        text=text,
                        x=x,
                        y=y,
                        z=z,
                        showarrow=False,
                        yshift=plot_config["label_offset"],
                    )
                )

        if show_distances:
            for (i, j), length in self.bond_lengths.items():
                labels.append(
                    dict(
                        text=f"{length:.2f}",
                        x=(self.x[i] + self.x[j]) / 2,
                        y=(self.y[i] + self.y[j]) / 2,
                        z=(self.z[i] + self.z[j]) / 2,
                        showarrow=False,
                        yshift=plot_config["label_offset"],
                        font=dict(color=plot_config["bond_label_color"]),
                    )
                )

        return labels

    def _create_menu_buttons(self, plot_config: PlotConfig) -> List[dict]:
        """Creates menu buttons for toggling different label displays."""
        return [
            dict(
                label="Elements",
                method="relayout",
                args=[
                    {
                        "scene.annotations": self._create_annotations(
                            plot_config, atom_label_type="element"
                        )
                    }
                ],
            ),
            dict(
                label="Elements & Distances",
                method="relayout",
                args=[
                    {
                        "scene.annotations": self._create_annotations(
                            plot_config,
                            atom_label_type="element",
                            show_distances=True,
                        )
                    }
                ],
            ),
            dict(
                label="IDs",
                method="relayout",
                args=[
                    {
                        "scene.annotations": self._create_annotations(
                            plot_config, atom_label_type="index"
                        )
                    }
                ],
            ),
            dict(
                label="IDs & Distances",
                method="relayout",
                args=[
                    {
                        "scene.annotations": self._create_annotations(
                            plot_config, atom_label_type="index", show_distances=True
                        )
                    }
                ],
            ),
            dict(
                label="Distances",
                method="relayout",
                args=[
                    {
                        "scene.annotations": self._create_annotations(
                            plot_config, show_distances=True
                        )
                    }
                ],
            ),
            dict(label="Clear", method="relayout", args=[{"scene.annotations": []}]),
        ]

    def to_plotly(
        self,
        config: Optional[PlotConfig] = None,
        title: Optional[str] = None,
    ) -> go.Figure:
        """Creates an interactive 3D visualization of the molecular structure.

        Generates a Plotly figure displaying atoms as spheres and bonds as lines, with
        configurable visual properties and interactive features for exploring the molecule's
        geometry.

        Args:
            config (Optional[PlotConfig]): Dictionary of visualization settings:

                - atom_size (int): Size of atom markers (default: 7)
                - atom_opacity (float): Opacity of atoms from 0.0-1.0 (default: 0.8)
                - atom_border_color (str): Color of atom borders (default: "lightgray")
                - atom_border_width (int): Width of atom borders (default: 2)
                - bond_color (str): Color of bonds between atoms (default: "grey")
                - bond_width (int): Width of bond lines (default: 2)
                - show_grid (bool): Whether to display axis grid (default: False)
                - label_offset (int): Vertical offset for labels (default: 15)
                - bond_label_color (str): Color of bond length labels (default: "steelblue")

            title (Optional[str]): Title displayed above visualization (default: None)

        Returns:
            go.Figure: Plotly figure object containing:

                - 3D scatter plot of atoms colored by element
                - Line segments representing bonds
                - Interactive menu for toggling labels:

                    - Element symbols
                    - Atom indices
                    - Bond lengths
                    - Combinations of the above

                - Mouse controls for rotation, zoom, and pan

        Examples:
            >>> fig = mg.to_plotly(config={
            ...     "atom_size": 12,
            ...     "atom_opacity": 0.9,
            ...     "bond_color": "black",
            ...     "show_grid": True,
            ...     "label_offset": 20
            ... })
            >>> fig.show()
        """
        logger.debug("Creating Plotly figure")

        if not self.elements:
            logger.warning("No molecular data to visualize")
            return go.Figure()

        plot_config = DEFAULT_PLOT_CONFIG.copy()
        if config:
            plot_config.update(config)

        # Create figure
        data = [
            self._create_atom_trace(plot_config),
            self._create_bond_trace(plot_config),
        ]

        # Configure axis parameters
        axis_params = dict(
            showgrid=plot_config["show_grid"],
            showbackground=False,
            showticklabels=False,
            zeroline=False,
            title={"text": ""},  # Hide axis titles
            showspikes=False,
        )

        # Build layout
        layout = dict(
            title=title,
            scene=dict(
                xaxis=axis_params,
                yaxis=axis_params,
                zaxis=axis_params,
                annotations=self._create_annotations(plot_config, atom_label_type="element"),
            ),
            margin=dict(r=0, l=0, b=0, t=0 if not title else 40),
            showlegend=False,
            updatemenus=[
                dict(
                    buttons=self._create_menu_buttons(plot_config),
                    direction="down",
                    xanchor="left",
                    yanchor="top",
                    pad=dict(r=10, t=10),
                )
            ],
        )

        logger.debug("Plotly figure created successfully")
        return go.Figure(data=data, layout=layout)

    def to_networkx(self) -> nx.Graph:
        """Convert the current graph representation to a NetworkX graph.

        Creates a NetworkX graph from the stored adjacency list, adding node and edge attributes.
        Node attributes include the chemical element and xyz coordinates, while edge attributes
        store bond lengths.

        Returns:
            nx.Graph: NetworkX graph containing molecular structure with node attributes 'element'
                and 'xyz' (coordinates), and edge attribute 'length' (bond length).

        Examples:
            >>> mg = MolGraph()
            >>> mg.read_xyz('molecule.xyz')
            >>> G = mg.to_networkx()
            >>> print(G.nodes[0]['element'])
            'C'
            >>> print(G.edges[(0, 1)]['length'])
            1.54
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
        edge_attrs = {edge: {"length": length} for edge, length in self.bond_lengths.items()}
        nx.set_edge_attributes(G, edge_attrs)
        logger.debug("NetworkX graph created successfully")

        return G

    def _generate_adjacency_list(self) -> None:
        """Generates the adjacency list and matrix based on atomic positions and radii.

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

            self.adj_matrix = np.logical_and(0.1 < distances, distance_bond > distances).astype(int)

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
        """Create an iterator yielding all bonds in the molecular graph.

        Iterates through the adjacency list to generate all unique pairs of bonded atoms
        in the structure, ensuring each bond is only yielded once.

        Yields:
            tuple[int, int]: A pair of atom indices representing a bonded atomic pair,
                where each index corresponds to the atom's position in the molecular structure.

        Examples:
            >>> mg = MolGraph()
            >>> mg.read_xyz('molecule.xyz')
            >>> for atom1, atom2 in mg.edges():
            ...     print(f"Bond between atoms {atom1} and {atom2}")
            Bond between atoms 0 and 1
            Bond between atoms 1 and 2
        """
        edges = set()
        for node, neighbours in self.adj_list.items():
            for neighbour in neighbours:
                bond = frozenset([node, neighbour])
                if bond in edges:
                    continue
                edges.add(bond)
                yield node, neighbour

    def formula(self) -> str:
        """Return the molecular formula of the structure in Hill notation.

        Generates a string representation of the molecular formula following Hill
        notation conventions, where carbon and hydrogen (if present) are listed first,
        followed by other elements in alphabetical order. Numbers are added as
        subscripts for elements appearing multiple times.

        Returns:
            str: The molecular formula string with elements ordered by Hill notation
                (e.g., 'C2H6O' for ethanol).

        Examples:
            >>> mg = MolGraph()
            >>> mg.read_xyz('ethanol.xyz')
            >>> print(mg.formula())
            'C2H6O'
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
        """Returns the number of atoms in the molecule.

        Returns:
            int: Number of atoms
        """
        return len(self.elements)

    def __getitem__(self, position: int) -> Tuple[str, Tuple[float, float, float]]:
        """Returns the element and coordinates for the atom at the given position.

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
        """Creates an iterator over all atoms in the molecule.

        Returns:
            Iterator[Tuple[str, Tuple[float, float, float]]]: Iterator yielding
            tuples of (element, coordinates) for each atom
        """
        return (self[i] for i in range(len(self)))

    def __repr__(self) -> str:
        """Creates a simplified string representation of the molecular graph.

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
