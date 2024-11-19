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
    Literal,
    Optional,
    Set,
    Tuple,
    TypedDict,
)

import plotly.graph_objects as go

from xyz2graph.geometry import Point3D

from .molecule import Atom, Bond


if TYPE_CHECKING:
    from .graph import MolGraph


class ButtonAppearance(TypedDict):
    """Configuration for button appearance."""

    direction: str
    pad_r: int
    pad_t: int
    showactive: bool
    xanchor: str
    yanchor: str


class ButtonPosition(TypedDict):
    """Position and title for a button group."""

    button_y: float
    label_y: float
    title: str


class ButtonsCommonConfig(TypedDict):
    """Common settings for all buttons."""

    button_x: float
    label_x: float


class ButtonsLayoutConfig(TypedDict):
    """Complete button configuration."""

    appearance: ButtonAppearance  # Common appearance settings
    traces: Dict[Literal["atoms", "bonds"], ButtonPosition]
    labels: Dict[Literal["atoms", "bonds"], ButtonPosition]
    background: ButtonPosition
    common: ButtonsCommonConfig


class StyleConfig(TypedDict, total=False):
    """Visual style configuration."""

    atom_size: int
    atom_opacity: float
    atom_border_color: str
    atom_border_width: int
    bond_color: str
    bond_width: int
    label_offset: float
    background_color: str
    axis_visible: bool
    title: str


class SceneConfig(TypedDict):
    """3D scene configuration."""

    showbackground: bool
    showticklabels: bool
    zeroline: bool
    showspikes: bool
    title: str


class VisualizationConfig(TypedDict):
    """Complete visualization configuration."""

    buttons: ButtonsLayoutConfig  # Button positions and labels
    style: StyleConfig  # Visual appearance settings
    scene: SceneConfig  # 3D scene settings


# Default configurations
DEFAULT_BUTTON_APPEARANCE = ButtonAppearance(
    direction="right", pad_r=10, pad_t=10, showactive=True, xanchor="left", yanchor="top"
)

DEFAULT_BUTTONS = ButtonsLayoutConfig(
    appearance=DEFAULT_BUTTON_APPEARANCE,
    traces={
        "atoms": {"button_y": 1.175, "label_y": 1.155, "title": "Atoms"},
        "bonds": {"button_y": 1.065, "label_y": 1.045, "title": "Bonds"},
    },
    labels={
        "atoms": {"button_y": 1.12, "label_y": 1.1, "title": "Atom Labels"},
        "bonds": {"button_y": 1.01, "label_y": 0.99, "title": "Bond Lengths"},
    },
    background={"button_y": 0.955, "label_y": 0.935, "title": "Background"},
    common={"button_x": 0.07, "label_x": 0.0},
)

DEFAULT_STYLE = StyleConfig(
    atom_size=7,
    atom_opacity=0.8,
    atom_border_color="lightgray",
    atom_border_width=2,
    bond_color="grey",
    bond_width=2,
    label_offset=0.3,
    background_color="white",
    axis_visible=True,
    title="",
)

DEFAULT_SCENE = SceneConfig(
    showbackground=True, showticklabels=False, zeroline=False, showspikes=False, title=""
)

DEFAULT_CONFIG = VisualizationConfig(
    buttons=DEFAULT_BUTTONS, style=DEFAULT_STYLE, scene=DEFAULT_SCENE
)


class ConfigurationManager:
    """Manages visualization configuration."""

    def __init__(self, config: Optional[VisualizationConfig] = None) -> None:
        """Initialize configuration manager.

        Args:
            config: Optional custom configuration. Uses DEFAULT_CONFIG if None.
        """
        self.config = config or DEFAULT_CONFIG

    def get_buttons_layout(self) -> ButtonsLayoutConfig:
        """Get button layout configuration.

        Returns:
            Button layout configuration dictionary.
        """
        return self.config["buttons"]

    def get_style(self) -> StyleConfig:
        """Get style configuration.

        Returns:
            Style configuration dictionary.
        """
        return self.config["style"]

    def get_scene(self) -> SceneConfig:
        """Get scene configuration.

        Returns:
            Scene configuration dictionary.
        """
        return self.config["scene"]

    @property
    def atom_style(self) -> Dict[str, Any]:
        """Get atom-specific style settings.

        Returns:
            Dictionary of atom style parameters.
        """
        style = self.get_style()
        return {
            "size": style["atom_size"],
            "opacity": style["atom_opacity"],
            "border_color": style["atom_border_color"],
            "border_width": style["atom_border_width"],
        }

    @property
    def bond_style(self) -> Dict[str, Any]:
        """Get bond-specific style settings.

        Returns:
            Dictionary of bond style parameters.
        """
        style = self.get_style()
        return {"color": style["bond_color"], "width": style["bond_width"]}

    @property
    def scene_params(self) -> Dict[str, Any]:
        """Get scene parameters for Plotly.

        Returns:
            Dictionary of scene parameters.
        """
        scene = self.get_scene()
        return {
            "xaxis": scene,
            "yaxis": scene,
            "zaxis": scene,
        }


class LabelType(Enum):
    """Types of labels in molecular visualization."""

    ATOM_INDEX = auto()
    ATOM_ELEMENT = auto()
    BOND_LENGTH = auto()


@dataclass
class BaseLabel(ABC):
    """Base class for molecular labels."""

    def __post_init__(self) -> None:
        """Initialize label attributes after instance creation."""
        self.show = False
        self.offset = 15
        self.color = "black"

    @abstractmethod
    def get_position(self) -> Point3D:
        """Get label position.

        Returns:
            3D point coordinates for label placement.
        """
        pass

    @abstractmethod
    def get_text(self) -> str:
        """Get label text.

        Returns:
            Text to display in the label.
        """
        pass

    def to_annotation(self) -> dict:
        """Convert to Plotly annotation.

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
            "yshift": self.offset,
            "font": {"size": 12, "color": self.color},
            "visible": self.show,
        }

    @abstractmethod
    def get_type(self) -> LabelType:
        """Get label type.

        Returns:
            Enumerated label type.
        """
        pass


@dataclass
class AtomLabel(BaseLabel):
    """Label for atoms."""

    atom: Atom
    label_type: LabelType = field(default=LabelType.ATOM_ELEMENT)

    def __post_init__(self) -> None:
        """Initialize atom label attributes."""
        super().__post_init__()

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

    def __post_init__(self) -> None:
        """Initialize bond label attributes."""
        super().__post_init__()

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
    config: VisualizationConfig
    selected_atoms: Set[int] = field(default_factory=set)

    def to_trace(self) -> go.Scatter3d:
        """Convert atoms to Plotly trace.

        Returns:
            Scatter3d trace representing atoms.
        """
        x = [atom.x for atom in self.atoms]
        y = [atom.y for atom in self.atoms]
        z = [atom.z for atom in self.atoms]

        colors = [self.colors.get(atom.element, "pink") for atom in self.atoms]
        sizes = [
            self.config["style"]["atom_size"] * (2 if atom.index in self.selected_atoms else 1)
            for atom in self.atoms
        ]

        hover_text = [
            f"{atom.index}. {atom.element} ({atom.x:.2f}, {atom.y:.2f}, {atom.z:.2f})"
            for atom in self.atoms
        ]

        return go.Scatter3d(
            x=x,
            y=y,
            z=z,
            mode="markers",
            marker=dict(
                size=sizes,
                color=colors,
                opacity=self.config["style"]["atom_opacity"],
                line=dict(
                    color=self.config["style"]["atom_border_color"],
                    width=self.config["style"]["atom_border_width"],
                ),
            ),
            hovertext=hover_text,
            hoverinfo="text",
            name="atoms",
        )


@dataclass
class BondComponent(VisualizationComponent):
    """Component for visualizing bonds."""

    bonds: List[Bond]
    config: VisualizationConfig
    selected_bonds: Set[int] = field(default_factory=set)

    def to_trace(self) -> go.Scatter3d:
        """Convert bonds to Plotly trace.

        Returns:
            Scatter3d trace representing bonds.
        """
        x, y, z = [], [], []

        for bond in self.bonds:
            x.extend([bond.atom1.x, bond.atom2.x, None])
            y.extend([bond.atom1.y, bond.atom2.y, None])
            z.extend([bond.atom1.z, bond.atom2.z, None])

        return go.Scatter3d(
            x=x,
            y=y,
            z=z,
            mode="lines",
            line=dict(
                color=self.config["style"]["bond_color"], width=self.config["style"]["bond_width"]
            ),
            hoverinfo="none",
            name="bonds",
        )


class ButtonFactory:
    """Factory for creating visualization control buttons."""

    def __init__(self, buttons_config: ButtonsLayoutConfig) -> None:
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

    def create_atom_trace_buttons(self) -> List[Dict[str, Any]]:
        """Create buttons for controlling atom trace visibility.

        Returns:
            List of button configurations for controlling atom trace visibility.
        """
        return [
            {
                "label": "Show",
                "method": "restyle",
                "args": [{"visible": True}, [0]],
            },
            {
                "label": "Hide",
                "method": "restyle",
                "args": [{"visible": False}, [0]],
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

    def create_bond_trace_buttons(self) -> List[Dict[str, Any]]:
        """Create buttons for controlling bond trace visibility.

        Returns:
            List of button configurations for controlling bond trace visibility.
        """
        return [
            {
                "label": "Show",
                "method": "restyle",
                "args": [{"visible": True}, [1]],
            },
            {
                "label": "Hide",
                "method": "restyle",
                "args": [{"visible": False}, [1]],
            },
        ]

    def create_background_buttons(self) -> List[Dict[str, Any]]:
        """Create buttons for controlling axis background visibility.

        Returns:
            List of button configurations for controlling background visibility.
        """
        return [
            {
                "label": "Show",
                "method": "relayout",
                "args": [
                    {
                        "scene.xaxis.showbackground": True,
                        "scene.yaxis.showbackground": True,
                        "scene.zaxis.showbackground": True,
                    }
                ],
            },
            {
                "label": "Hide",
                "method": "relayout",
                "args": [
                    {
                        "scene.xaxis.showbackground": False,
                        "scene.yaxis.showbackground": False,
                        "scene.zaxis.showbackground": False,
                    }
                ],
            },
        ]

    def create_toggle_button_group(
        self,
        position: ButtonPosition,
        buttons_config: List[Dict[str, Any]],
        active_index: int = 0,
    ) -> Tuple[Dict[str, Any], Dict[str, Any]]:
        """Create a button group with title annotation.

        Args:
            position: Position and title configuration for the button group.
            buttons_config: List of button configurations.
            active_index: Index of initially active button.

        Returns:
            Tuple containing button group and title annotation configurations.
        """
        appearance = self.config["appearance"]
        common = self.config["common"]

        button_group = {
            "buttons": buttons_config,
            "direction": appearance["direction"],
            "pad": {"r": appearance["pad_r"], "t": appearance["pad_t"]},
            "showactive": appearance["showactive"],
            "active": active_index,
            "x": common["button_x"],
            "y": position["button_y"],
            "xanchor": appearance["xanchor"],
            "yanchor": appearance["yanchor"],
            "type": "buttons",
        }

        title_annotation = {
            "text": position["title"],
            "x": common["label_x"],
            "xref": "paper",
            "y": position["label_y"],
            "yref": "paper",
            "align": "left",
            "showarrow": False,
        }

        return button_group, title_annotation


class VisualizationControls:
    """Manages visualization control buttons and their interactions."""

    def __init__(self, n_atoms: int, n_bonds: int, buttons_config: ButtonsLayoutConfig) -> None:
        """Initialize visualization controls.

        Args:
            n_atoms: Number of atoms in the molecule.
            n_bonds: Number of bonds in the molecule.
            buttons_config: Configuration for button layout and appearance.
        """
        self.n_atoms = n_atoms
        self.n_bonds = n_bonds
        self.config = buttons_config
        self.button_factory = ButtonFactory(buttons_config)
        self.element_indices: List[int] = []
        self.index_indices: List[int] = []
        self.bond_indices: List[int] = []

    def update_annotation_indices(
        self, element_indices: List[int], index_indices: List[int], bond_indices: List[int]
    ) -> None:
        """Update the annotation indices for proper button targeting.

        Args:
            element_indices: List of indices for element labels.
            index_indices: List of indices for atom index labels.
            bond_indices: List of indices for bond labels.
        """
        self.element_indices = element_indices
        self.index_indices = index_indices
        self.bond_indices = bond_indices

    def create_all_controls(self) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
        """Create all control buttons and their annotations.

        Returns:
            Tuple containing lists of button configurations and annotations.
        """
        buttons = []
        annotations = []

        # Atom trace controls
        atom_trace_group, atom_trace_annotation = self.button_factory.create_toggle_button_group(
            position=self.config["traces"]["atoms"],
            buttons_config=self.button_factory.create_atom_trace_buttons(),
        )
        buttons.append(atom_trace_group)
        annotations.append(atom_trace_annotation)

        # Atom label controls
        if self.n_atoms > 0:
            atom_label_group, atom_label_annotation = (
                self.button_factory.create_toggle_button_group(
                    position=self.config["labels"]["atoms"],
                    buttons_config=self.button_factory.create_atom_label_buttons(
                        self.element_indices, self.index_indices
                    ),
                    active_index=2,  # Initially hidden
                )
            )
            buttons.append(atom_label_group)
            annotations.append(atom_label_annotation)

        # Bond trace controls
        if self.n_bonds > 0:
            bond_trace_group, bond_trace_annotation = (
                self.button_factory.create_toggle_button_group(
                    position=self.config["traces"]["bonds"],
                    buttons_config=self.button_factory.create_bond_trace_buttons(),
                )
            )
            buttons.append(bond_trace_group)
            annotations.append(bond_trace_annotation)

            bond_label_group, bond_label_annotation = (
                self.button_factory.create_toggle_button_group(
                    position=self.config["labels"]["bonds"],
                    buttons_config=self.button_factory.create_bond_label_buttons(self.bond_indices),
                    active_index=1,  # Initially hidden
                )
            )
            buttons.append(bond_label_group)
            annotations.append(bond_label_annotation)

        # Background controls
        bg_group, bg_annotation = self.button_factory.create_toggle_button_group(
            position=self.config["background"],
            buttons_config=self.button_factory.create_background_buttons(),
        )
        buttons.append(bg_group)
        annotations.append(bg_annotation)

        return buttons, annotations


class VisualizationManager:
    """Manages the complete visualization of molecular structures."""

    def __init__(self, mol_graph: "MolGraph", config: Optional[VisualizationConfig] = None) -> None:
        """Initialize visualization manager.

        Args:
            mol_graph: Molecular graph to visualize.
            config: Optional custom visualization configuration.
        """
        self.mol_graph = mol_graph
        self.config = config or DEFAULT_CONFIG

        self.atom_component = AtomComponent(
            atoms=mol_graph.atoms,
            colors=mol_graph.cpk_colors,
            config=self.config,
            selected_atoms=set(),
        )

        self.bond_component = BondComponent(
            bonds=mol_graph.bonds,
            config=self.config,
            selected_bonds=set(),
        )

        self.atom_element_labels = [
            AtomLabel(atom=atom, label_type=LabelType.ATOM_ELEMENT) for atom in mol_graph.atoms
        ]

        self.atom_index_labels = [
            AtomLabel(atom=atom, label_type=LabelType.ATOM_INDEX) for atom in mol_graph.atoms
        ]

        self.bond_labels = [BondLabel(bond=bond) for bond in mol_graph.bonds]

        self.controls = VisualizationControls(
            n_atoms=len(mol_graph.atoms),
            n_bonds=len(mol_graph.bonds),
            buttons_config=self.config["buttons"],
        )

    def create_figure(self) -> go.Figure:
        """Create a complete interactive figure for molecular visualization.

        Returns:
            Plotly Figure object containing the molecular visualization.
        """
        fig = go.Figure()

        # Add main traces
        fig.add_trace(self.atom_component.to_trace())
        fig.add_trace(self.bond_component.to_trace())

        # Prepare annotations
        element_annotations = [label.to_annotation() for label in self.atom_element_labels]
        index_annotations = [label.to_annotation() for label in self.atom_index_labels]
        bond_annotations = [label.to_annotation() for label in self.bond_labels]

        for ann in element_annotations + index_annotations + bond_annotations:
            ann["visible"] = False

        # Calculate annotation indices
        self.element_annotation_indices = list(range(len(element_annotations)))
        self.index_annotation_indices = list(
            range(len(element_annotations), len(element_annotations) + len(index_annotations))
        )
        self.bond_annotation_indices = list(
            range(
                len(element_annotations) + len(index_annotations),
                len(element_annotations) + len(index_annotations) + len(bond_annotations),
            )
        )

        self.controls.update_annotation_indices(
            element_indices=self.element_annotation_indices,
            index_indices=self.index_annotation_indices,
            bond_indices=self.bond_annotation_indices,
        )

        buttons, control_annotations = self.controls.create_all_controls()

        # Create watermark
        watermark = dict(
            text=(
                'created with <a href="https://zotko.github.io/xyz2graph/" '
                'style="color:#404040;">xyz2graph</a>'
            ),
            xref="paper",
            yref="paper",
            x=0.99,
            y=0.01,
            showarrow=False,
            font=dict(size=10, color="gray"),
            xanchor="right",
            yanchor="bottom",
        )

        # Set up axis parameters
        axis_params = {
            "showbackground": self.config["scene"]["showbackground"],
            "showticklabels": self.config["scene"]["showticklabels"],
            "zeroline": self.config["scene"]["zeroline"],
            "showspikes": self.config["scene"]["showspikes"],
            "title": {"text": self.config["scene"]["title"]},
        }

        # Update layout
        fig.update_layout(
            title=self.config["style"].get("title", ""),
            scene=dict(
                xaxis=axis_params,
                yaxis=axis_params,
                zaxis=axis_params,
                annotations=element_annotations + index_annotations + bond_annotations,
                camera=dict(
                    up=dict(x=0, y=0, z=1),
                    center=dict(x=0, y=0, z=0),
                    eye=dict(x=1.5, y=1.5, z=1.5),
                ),
            ),
            margin=dict(r=0, l=0, b=0, t=0 if not self.config["style"].get("title") else 40),
            showlegend=False,
            updatemenus=buttons,
            annotations=control_annotations + [watermark],
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
