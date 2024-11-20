"""Visualization module for molecular structures using Plotly.

This module provides classes and utilities for creating interactive 3D
visualizations of molecular structures. It includes support for atoms,
bonds, labels, and interactive controls.
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from enum import Enum, auto
from typing import (
    TYPE_CHECKING,
    Any,
    Dict,
    List,
    Mapping,
    Optional,
    Tuple,
    TypedDict,
)

import plotly.graph_objects as go

from .geometry import Point3D
from .molecule import Atom, Bond


if TYPE_CHECKING:
    from .graph import MolGraph


# Button-related configurations
class ButtonAppearance(TypedDict):
    """Common appearance settings for all buttons."""

    direction: str
    pad_r: int
    pad_t: int
    showactive: bool
    xanchor: str
    yanchor: str
    font_size: int  # Add this as it was hardcoded (10)


class ButtonPosition(TypedDict):
    """Position and title for a button group."""

    button_y: float
    label_y: float
    title: str


class ButtonsConfig(TypedDict):
    """Complete button configuration."""

    appearance: ButtonAppearance
    position: Dict[str, float]  # Common x positions
    groups: Mapping[str, ButtonPosition]  # Combines traces, labels, and background


# Component configurations
class LabelConfig(TypedDict):
    """Label display configuration."""

    offset: float
    font_size: int
    color: str
    visible: bool


class AtomConfig(TypedDict):
    """Atom visualization configuration."""

    size: int
    opacity: float
    border_color: str
    border_width: int
    hover_font_size: int


class BondConfig(TypedDict):
    """Bond visualization configuration."""

    color: str
    width: int


class LegendConfig(TypedDict):
    """Legend display configuration."""

    x: float
    y: float
    xanchor: str
    yanchor: str


# Scene configuration
class SceneConfig(TypedDict):
    """3D scene configuration."""

    showbackground: bool
    showticklabels: bool
    zeroline: bool
    showspikes: bool
    background_color: str  # Moved from StyleConfig
    axis_visible: bool  # Moved from StyleConfig
    legend: LegendConfig


# Main configuration structures
class StyleConfig(TypedDict):
    """Visual style configuration."""

    atoms: AtomConfig
    bonds: BondConfig
    labels: LabelConfig
    title: str


class VisualizationConfig(TypedDict):
    """Complete visualization configuration."""

    style: StyleConfig
    scene: SceneConfig
    buttons: ButtonsConfig


# Default configurations
DEFAULT_BUTTON_APPEARANCE = ButtonAppearance(
    direction="right",
    pad_r=5,
    pad_t=5,
    showactive=True,
    xanchor="left",
    yanchor="top",
    font_size=10,
)

DEFAULT_BUTTON_POSITION = {"button_x": 0.05, "label_x": 0.0}

DEFAULT_BUTTON_GROUPS: Dict[str, ButtonPosition] = {
    "atoms": {"button_y": 1.0, "label_y": 0.99, "title": "Atoms"},
    "atom_labels": {"button_y": 0.95, "label_y": 0.94, "title": "Atom Labels"},
    "bonds": {"button_y": 0.90, "label_y": 0.89, "title": "Bonds"},
    "bond_labels": {"button_y": 0.85, "label_y": 0.84, "title": "Bond Labels"},
    "background": {"button_y": 0.80, "label_y": 0.79, "title": "Background"},
}

DEFAULT_BUTTONS = ButtonsConfig(
    appearance=DEFAULT_BUTTON_APPEARANCE,
    position=DEFAULT_BUTTON_POSITION,
    groups=DEFAULT_BUTTON_GROUPS,
)

DEFAULT_LABEL_CONFIG = LabelConfig(offset=15, font_size=12, color="black", visible=False)

DEFAULT_ATOM_CONFIG = AtomConfig(
    size=7, opacity=0.8, border_color="lightgray", border_width=2, hover_font_size=12
)

DEFAULT_BOND_CONFIG = BondConfig(color="grey", width=2)

DEFAULT_SCENE_CONFIG = SceneConfig(
    showbackground=False,
    showticklabels=False,
    zeroline=False,
    showspikes=False,
    background_color="white",
    axis_visible=True,
    legend=LegendConfig(x=0.99, y=0.96, xanchor="right", yanchor="top"),
)

DEFAULT_STYLE = StyleConfig(
    atoms=DEFAULT_ATOM_CONFIG, bonds=DEFAULT_BOND_CONFIG, labels=DEFAULT_LABEL_CONFIG, title=""
)

DEFAULT_CONFIG = VisualizationConfig(
    style=DEFAULT_STYLE, scene=DEFAULT_SCENE_CONFIG, buttons=DEFAULT_BUTTONS
)

WATERMARK = {
    "text": (
        'created with <a href="https://zotko.github.io/xyz2graph/"'
        ' style="color:#404040;">xyz2graph</a>'
    ),
    "xref": "paper",
    "yref": "paper",
    "x": 0.99,
    "y": 0.01,
    "showarrow": False,
    "font": {"size": 10, "color": "gray"},
    "xanchor": "right",
    "yanchor": "bottom",
}


class ConfigurationManager:
    """Provides access to different configuration aspects."""

    def __init__(self, config: Optional[VisualizationConfig] = None) -> None:
        """Initialize configuration manager.

        Args:
            config: Custom visualization configuration. Uses DEFAULT_CONFIG if None.
        """
        self.config = config or DEFAULT_CONFIG

    @property
    def style(self) -> StyleConfig:
        """Get style configuration."""
        return self.config["style"]

    @property
    def scene(self) -> SceneConfig:
        """Get scene configuration."""
        return self.config["scene"]

    @property
    def buttons(self) -> ButtonsConfig:
        """Get buttons configuration."""
        return self.config["buttons"]

    @property
    def atoms(self) -> AtomConfig:
        """Get atom visualization configuration."""
        return self.style["atoms"]

    @property
    def bonds(self) -> BondConfig:
        """Get bond visualization configuration."""
        return self.style["bonds"]

    @property
    def labels(self) -> LabelConfig:
        """Get label display configuration."""
        return self.style["labels"]


class LabelType(Enum):
    """Types of labels in molecular visualization."""

    ATOM_INDEX = auto()
    ATOM_ELEMENT = auto()
    BOND_LENGTH = auto()


@dataclass
class BaseLabel(ABC):
    """Base class for molecular labels."""

    def to_annotation(self, config: LabelConfig = DEFAULT_LABEL_CONFIG) -> dict:
        """Convert label to Plotly annotation format.

        Args:
            config: Label display configuration.

        Returns:
            Dictionary containing Plotly annotation parameters.
        """
        pos = self.get_position()

        return {
            "text": self.get_text(),
            "x": pos.x,
            "y": pos.y,
            "z": pos.z,
            "showarrow": False,
            "yshift": config["offset"],
            "font": {"size": config["font_size"], "color": config["color"]},
            "visible": config["visible"],
        }

    @abstractmethod
    def get_position(self) -> Point3D:
        """Get label position in 3D space."""
        pass

    @abstractmethod
    def get_text(self) -> str:
        """Get label text content."""
        pass

    @abstractmethod
    def get_type(self) -> LabelType:
        """Get the type of this label."""
        pass


@dataclass
class AtomLabel(BaseLabel):
    """Label for atoms."""

    atom: Atom
    label_type: LabelType = field(default=LabelType.ATOM_ELEMENT)

    def get_position(self) -> Point3D:
        """Get atom position for label placement."""
        return self.atom.position

    def get_text(self) -> str:
        """Get atom label text based on label type."""
        if self.label_type == LabelType.ATOM_ELEMENT:
            return self.atom.element
        return str(self.atom.index)

    def get_type(self) -> LabelType:
        """Get atom label type."""
        return self.label_type


@dataclass
class BondLabel(BaseLabel):
    """Label for bonds."""

    bond: Bond

    def get_position(self) -> Point3D:
        """Get bond midpoint position for label placement."""
        return Point3D.midpoint(self.bond.atom1.position, self.bond.atom2.position)

    def get_text(self) -> str:
        """Get bond length formatted as text."""
        return f"{self.bond.length:.2f}"

    def get_type(self) -> LabelType:
        """Get bond label type."""
        return LabelType.BOND_LENGTH


class VisualizationComponent(ABC):
    """Base class for visualization components."""

    @abstractmethod
    def to_trace(self) -> go.Scatter3d:
        """Convert component to Plotly trace.

        Returns:
            Plotly Scatter3d trace object.
        """
        pass


@dataclass
class AtomComponent(VisualizationComponent):
    """Component for visualizing atoms."""

    atoms: List[Atom]
    colors: Dict[str, str]

    def _create_hover_text(self) -> List[str]:
        """Create hover text for atoms."""
        return [
            f"{atom.index}. {atom.element} ({atom.x:.2f}, {atom.y:.2f}, {atom.z:.2f})"
            for atom in self.atoms
        ]

    def _get_coordinates(self) -> Tuple[List[float], List[float], List[float]]:
        """Get atom coordinates."""
        return (
            [atom.x for atom in self.atoms],
            [atom.y for atom in self.atoms],
            [atom.z for atom in self.atoms],
        )

    def _get_colors(self) -> List[str]:
        """Get atom colors."""
        return [self.colors.get(atom.element, "pink") for atom in self.atoms]

    def to_trace(self, config: AtomConfig = DEFAULT_ATOM_CONFIG) -> List[go.Scatter3d]:
        """Convert atoms to Plotly traces, grouped by element.

        Args:
            config: Atom visualization configuration.

        Returns:
            List of Scatter3d traces, one per unique element type.
        """
        # Group atoms by element
        element_groups: Dict[str, List[Atom]] = {}
        for atom in self.atoms:
            if atom.element not in element_groups:
                element_groups[atom.element] = []
            element_groups[atom.element].append(atom)

        traces = []
        for element, atoms in element_groups.items():
            # Get coordinates for this element group
            x = [atom.x for atom in atoms]
            y = [atom.y for atom in atoms]
            z = [atom.z for atom in atoms]

            # Create hover text for this group
            hover_text = [
                f"{atom.index}. {atom.element} ({atom.x:.2f}, {atom.y:.2f}, {atom.z:.2f})"
                for atom in atoms
            ]

            # Create trace for this element
            trace = go.Scatter3d(
                x=x,
                y=y,
                z=z,
                mode="markers",
                marker=dict(
                    size=config["size"],
                    color=self.colors.get(element, "pink"),
                    opacity=config["opacity"],
                    line=dict(
                        color=config["border_color"],
                        width=config["border_width"],
                    ),
                ),
                hovertext=hover_text,
                hoverinfo="text",
                name=element,
            )
            traces.append(trace)

        return traces


@dataclass
class BondComponent(VisualizationComponent):
    """Component for visualizing bonds."""

    bonds: List[Bond]

    def _get_coordinates(
        self,
    ) -> Tuple[List[Optional[float]], List[Optional[float]], List[Optional[float]]]:
        """Get bond coordinates including None separators."""
        x, y, z = [], [], []

        for bond in self.bonds:
            x.extend([bond.atom1.x, bond.atom2.x, None])
            y.extend([bond.atom1.y, bond.atom2.y, None])
            z.extend([bond.atom1.z, bond.atom2.z, None])

        return x, y, z

    def to_trace(self, config: BondConfig = DEFAULT_BOND_CONFIG) -> go.Scatter3d:
        """Convert bonds to Plotly trace.

        Args:
            config: Bond visualization configuration.

        Returns:
            Scatter3d trace representing bonds.
        """
        x, y, z = self._get_coordinates()

        return go.Scatter3d(
            x=x,
            y=y,
            z=z,
            mode="lines",
            line=dict(color=config["color"], width=config["width"]),
            hoverinfo="none",
            name="bonds",
        )


class ButtonFactory:
    """Factory for creating visualization control buttons."""

    def __init__(self, buttons_config: ButtonsConfig) -> None:
        """Initialize button factory.

        Args:
            buttons_config: Configuration for button layout and appearance.
        """
        self.config = buttons_config

    def create_atom_label_buttons(
        self, element_indices: List[int], index_indices: List[int]
    ) -> List[Dict[str, Any]]:
        """Create buttons for controlling atom label visibility.

        Args:
            element_indices: Indices of element labels.
            index_indices: Indices of atom index labels.

        Returns:
            List of button configurations for controlling atom label visibility.
        """
        return [
            {
                "label": "ID",
                "method": "update",
                "args": [
                    {},
                    {
                        **{f"scene.annotations[{i}].visible": True for i in index_indices},
                        **{f"scene.annotations[{i}].visible": False for i in element_indices},
                    },
                ],
            },
            {
                "label": "Element",
                "method": "update",
                "args": [
                    {},
                    {
                        **{f"scene.annotations[{i}].visible": False for i in index_indices},
                        **{f"scene.annotations[{i}].visible": True for i in element_indices},
                    },
                ],
            },
            {
                "label": "Hide",
                "method": "update",
                "args": [
                    {},
                    {
                        **{f"scene.annotations[{i}].visible": False for i in index_indices},
                        **{f"scene.annotations[{i}].visible": False for i in element_indices},
                    },
                ],
            },
        ]

    def create_bond_label_buttons(self, bond_indices: List[int]) -> List[Dict[str, Any]]:
        """Create buttons for controlling bond label visibility.

        Args:
            bond_indices: Indices of bond labels.

        Returns:
            List of button configurations for controlling bond label visibility.
        """
        return [
            {
                "label": "Show",
                "method": "update",
                "args": [{}, {**{f"scene.annotations[{i}].visible": True for i in bond_indices}}],
            },
            {
                "label": "Hide",
                "method": "update",
                "args": [{}, {**{f"scene.annotations[{i}].visible": False for i in bond_indices}}],
            },
        ]

    def create_atom_trace_buttons(self, atom_trace_indices: List[int]) -> List[Dict[str, Any]]:
        """Create buttons for controlling all atom traces visibility together.

        Args:
            atom_trace_indices: List of trace indices for all atom elements

        Returns:
            Button configurations for controlling all atoms visibility together
        """
        return [
            {
                "label": "Show",
                "method": "restyle",
                "args": [{"visible": True}, atom_trace_indices],
            },
            {
                "label": "Hide",
                "method": "restyle",
                "args": [{"visible": False}, atom_trace_indices],
            },
        ]

    def create_bond_trace_buttons(self, bond_trace_indices: List[int]) -> List[Dict[str, Any]]:
        """Create buttons for controlling bond trace visibility."""
        return [
            {
                "label": "Show",
                "method": "restyle",
                "args": [{"visible": True}, bond_trace_indices],
            },
            {
                "label": "Hide",
                "method": "restyle",
                "args": [{"visible": False}, bond_trace_indices],
            },
        ]

    def create_background_buttons(self) -> List[Dict[str, Any]]:
        """Create buttons for controlling axis background visibility."""

        def make_args(show: bool) -> Dict[str, bool]:
            return {f"scene.{axis}.showbackground": show for axis in ["xaxis", "yaxis", "zaxis"]}

        return [
            {
                "label": "Show",
                "method": "relayout",
                "args": [make_args(True)],
            },
            {
                "label": "Hide",
                "method": "relayout",
                "args": [make_args(False)],
            },
        ]

    def create_button_group(
        self, group_name: str, buttons: List[Dict[str, Any]], active_index: int = 0
    ) -> Tuple[Dict[str, Any], Dict[str, Any]]:
        """Create a button group with title annotation.

        Args:
            group_name: Name of the button group
            buttons: List of button configurations
            active_index: Index of initially active button

        Returns:
            Tuple of (button group configuration, title annotation)
        """
        group_position = self.config["groups"][group_name]
        appearance = self.config["appearance"]
        position = self.config["position"]

        button_group = {
            "buttons": buttons,
            "direction": appearance["direction"],
            "pad": {"r": appearance["pad_r"], "t": appearance["pad_t"]},
            "showactive": appearance["showactive"],
            "active": active_index,
            "x": position["button_x"],
            "y": group_position["button_y"],
            "xanchor": appearance["xanchor"],
            "yanchor": appearance["yanchor"],
            "type": "buttons",
            "font": {"size": appearance["font_size"]},
        }

        title_annotation = {
            "text": group_position["title"],
            "x": position["label_x"],
            "xref": "paper",
            "y": group_position["label_y"],
            "yref": "paper",
            "align": "left",
            "showarrow": False,
            "font": {"size": appearance["font_size"]},
        }

        return button_group, title_annotation


class VisualizationControls:
    """Manages visualization control buttons and their interactions."""

    def __init__(self, n_atoms: int, n_bonds: int, buttons_config: ButtonsConfig) -> None:
        """Initialize visualization controls.

        Args:
            n_atoms: Number of atoms in the molecule.
            n_bonds: Number of bonds in the molecule.
            buttons_config: Configuration for buttons layout and appearance.
        """
        self.n_atoms = n_atoms
        self.n_bonds = n_bonds
        self.button_factory = ButtonFactory(buttons_config)

        # Indices for managing traces and annotations
        self.atom_trace_indices: List[int] = []
        self.bond_trace_indices: List[int] = []

        # Indices for managing annotations
        self.element_indices: List[int] = []
        self.index_indices: List[int] = []
        self.bond_indices: List[int] = []

    def update_annotation_indices(
        self, element_indices: List[int], index_indices: List[int], bond_indices: List[int]
    ) -> None:
        """Update the annotation indices for proper button targeting."""
        self.element_indices = element_indices
        self.index_indices = index_indices
        self.bond_indices = bond_indices

    def update_trace_indices(
        self, atom_trace_indices: List[int], bond_trace_indices: List[int]
    ) -> None:
        """Update the trace indices for proper button targeting."""
        self.atom_trace_indices = atom_trace_indices
        self.bond_trace_indices = bond_trace_indices

    def create_all_controls(self) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
        """Create all control buttons and their annotations.

        Returns:
            Tuple containing lists of button configurations and annotations.
        """
        buttons = []
        annotations = []

        # Atom trace controls
        atom_trace_buttons = self.button_factory.create_atom_trace_buttons(self.atom_trace_indices)
        button_group, annotation = self.button_factory.create_button_group(
            "atoms", atom_trace_buttons
        )
        buttons.append(button_group)
        annotations.append(annotation)

        # Atom label controls
        if self.n_atoms > 0:
            atom_label_buttons = self.button_factory.create_atom_label_buttons(
                self.element_indices, self.index_indices
            )
            button_group, annotation = self.button_factory.create_button_group(
                "atom_labels",
                atom_label_buttons,
                active_index=2,  # Initially hidden
            )
            buttons.append(button_group)
            annotations.append(annotation)

        # Bond trace controls
        if self.n_bonds > 0:
            bond_trace_buttons = self.button_factory.create_bond_trace_buttons(
                self.bond_trace_indices
            )
            button_group, annotation = self.button_factory.create_button_group(
                "bonds", bond_trace_buttons
            )
            buttons.append(button_group)
            annotations.append(annotation)

            bond_label_buttons = self.button_factory.create_bond_label_buttons(self.bond_indices)
            button_group, annotation = self.button_factory.create_button_group(
                "bond_labels",
                bond_label_buttons,
                active_index=1,  # Initially hidden
            )
            buttons.append(button_group)
            annotations.append(annotation)

        # Background controls
        background_buttons = self.button_factory.create_background_buttons()
        button_group, annotation = self.button_factory.create_button_group(
            "background",
            background_buttons,
            active_index=1,  # Initially hidden
        )
        buttons.append(button_group)
        annotations.append(annotation)

        return buttons, annotations


class VisualizationManager:
    """Manages the complete visualization of molecular structures."""

    def __init__(self, mol_graph: "MolGraph", config: Optional[VisualizationConfig] = None) -> None:
        """Initialize visualization manager.

        Args:
            mol_graph: Molecular graph to visualize.
            config: Optional custom visualization configuration.
        """
        # Initialize configuration manager
        self.config_manager = ConfigurationManager(config)

        # Initialize visualization components
        self.atom_component = AtomComponent(atoms=mol_graph.atoms, colors=mol_graph.cpk_colors)

        self.bond_component = BondComponent(bonds=mol_graph.bonds)

        # Initialize labels
        self.atom_element_labels = [
            AtomLabel(atom=atom, label_type=LabelType.ATOM_ELEMENT) for atom in mol_graph.atoms
        ]

        self.atom_index_labels = [
            AtomLabel(atom=atom, label_type=LabelType.ATOM_INDEX) for atom in mol_graph.atoms
        ]

        self.bond_labels = [BondLabel(bond=bond) for bond in mol_graph.bonds]

        # Initialize controls
        self.controls = VisualizationControls(
            n_atoms=len(mol_graph.atoms),
            n_bonds=len(mol_graph.bonds),
            buttons_config=self.config_manager.buttons,
        )

    def _create_annotations(self) -> Tuple[List[Dict[str, Any]], List[int], List[int], List[int]]:
        """Create and configure all annotations.

        Returns:
            Tuple containing:
            - List of all annotations
            - List of element label indices
            - List of atom index label indices
            - List of bond label indices
        """
        # Create annotations with configuration
        element_annotations = [
            label.to_annotation(config=self.config_manager.labels)
            for label in self.atom_element_labels
        ]

        index_annotations = [
            label.to_annotation(config=self.config_manager.labels)
            for label in self.atom_index_labels
        ]

        bond_annotations = [
            label.to_annotation(config=self.config_manager.labels) for label in self.bond_labels
        ]

        # Set initial visibility to False
        for ann in element_annotations + index_annotations + bond_annotations:
            ann["visible"] = False

        # Calculate annotation indices
        element_indices = list(range(len(element_annotations)))
        index_indices = list(
            range(len(element_annotations), len(element_annotations) + len(index_annotations))
        )
        bond_indices = list(
            range(
                len(element_annotations) + len(index_annotations),
                len(element_annotations) + len(index_annotations) + len(bond_annotations),
            )
        )

        return (
            element_annotations + index_annotations + bond_annotations,
            element_indices,
            index_indices,
            bond_indices,
        )

    def _create_axis_config(self) -> Dict[str, Dict[str, Any]]:
        """Create axis configuration for all three axes.

        Returns:
            Dictionary containing axis configurations.
        """
        axis_config = {
            "showbackground": self.config_manager.scene["showbackground"],
            "showticklabels": self.config_manager.scene["showticklabels"],
            "zeroline": self.config_manager.scene["zeroline"],
            "showspikes": self.config_manager.scene["showspikes"],
            "visible": self.config_manager.scene["axis_visible"],
            "title": {"text": ""},  # Hide axis title
        }

        # Add background color if showbackground is True
        if self.config_manager.scene["showbackground"]:
            axis_config["backgroundcolor"] = self.config_manager.scene["background_color"]

        return {axis: axis_config.copy() for axis in ["xaxis", "yaxis", "zaxis"]}

    def create_figure(self) -> go.Figure:
        """Create a complete interactive figure for molecular visualization."""
        fig = go.Figure()

        # Add main traces with configuration
        atom_traces = self.atom_component.to_trace(config=self.config_manager.atoms)
        for trace in atom_traces:
            fig.add_trace(trace)

        bond_trace = self.bond_component.to_trace(config=self.config_manager.bonds)
        bond_trace.showlegend = False
        fig.add_trace(bond_trace)

        all_annotations, element_indices, index_indices, bond_indices = self._create_annotations()

        self.controls.update_annotation_indices(
            element_indices=element_indices, index_indices=index_indices, bond_indices=bond_indices
        )
        self.controls.update_trace_indices(
            atom_trace_indices=list(range(len(atom_traces))),
            bond_trace_indices=[len(atom_traces)],
        )

        buttons, control_annotations = self.controls.create_all_controls()
        axis_config = self._create_axis_config()
        top_margin = 40 if self.config_manager.style["title"] else 0

        fig.update_layout(
            title=self.config_manager.style["title"],
            scene={
                **axis_config,
                "annotations": all_annotations,
            },
            updatemenus=buttons,
            annotations=control_annotations + [WATERMARK],
            margin=dict(r=0, l=0, b=0, t=top_margin),
            legend=self.config_manager.scene["legend"],
        )

        return fig


def create_visualization(
    mol_graph: "MolGraph", config: Optional[VisualizationConfig] = None
) -> go.Figure:
    """Create a visualization of a molecular graph.

    Args:
        mol_graph: Molecular graph to visualize.
        config: Optional custom visualization configuration.

    Returns:
        Plotly Figure object containing the molecular visualization.
    """
    manager = VisualizationManager(mol_graph, config)
    return manager.create_figure()
