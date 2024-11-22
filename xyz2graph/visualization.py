"""Visualization module for molecular structures using Plotly.

This module provides classes and utilities for creating interactive 3D
visualizations of molecular structures. It includes support for atoms,
bonds, labels, and interactive controls.
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass
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

import numpy as np
import plotly.graph_objects as go

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
    font_size: int
    title_font_size: int


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
    itemsizing: str


# Scene configuration
class SceneConfig(TypedDict):
    """3D scene configuration."""

    showbackground: bool
    showticklabels: bool
    zeroline: bool
    showspikes: bool
    background_color: str
    axis_visible: bool
    legend: LegendConfig


# Main configuration structures
class StyleConfig(TypedDict):
    """Visual style configuration."""

    atoms: AtomConfig
    bonds: BondConfig
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
    title_font_size=12,
)

DEFAULT_BUTTON_POSITION = {"button_x": 0.04, "label_x": 0.0}

DEFAULT_BUTTON_GROUPS: Dict[str, ButtonPosition] = {
    "atoms": {"button_y": 1.0, "label_y": 0.99, "title": "Atoms"},
    "bonds": {"button_y": 0.95, "label_y": 0.94, "title": "Bonds"},
    "grid": {"button_y": 0.90, "label_y": 0.89, "title": "Grid"},
}

DEFAULT_BUTTONS = ButtonsConfig(
    appearance=DEFAULT_BUTTON_APPEARANCE,
    position=DEFAULT_BUTTON_POSITION,
    groups=DEFAULT_BUTTON_GROUPS,
)

DEFAULT_ATOM_CONFIG = AtomConfig(
    size=5, opacity=0.9, border_color="lightgray", border_width=2, hover_font_size=12
)

DEFAULT_BOND_CONFIG = BondConfig(color="grey", width=2)

DEFAULT_LEGEND_CONFIG = LegendConfig(
    x=0.99, y=0.96, xanchor="right", yanchor="top", itemsizing="constant"
)

DEFAULT_SCENE_CONFIG = SceneConfig(
    showbackground=True,
    showticklabels=False,
    zeroline=False,
    showspikes=False,
    background_color="white",
    axis_visible=True,
    legend=DEFAULT_LEGEND_CONFIG,
)

DEFAULT_STYLE = StyleConfig(atoms=DEFAULT_ATOM_CONFIG, bonds=DEFAULT_BOND_CONFIG, title="")

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
class AtomShape:
    """A class for managing atom shapes in the visualization.

    Handles the conversion of atomic data into Plotly visualization traces,
    including grouping atoms by element and applying visual styling.
    """

    atoms: List[Atom]
    colors: Dict[str, str]

    def to_trace(self, config: AtomConfig = DEFAULT_ATOM_CONFIG) -> List[go.Scatter3d]:
        """Convert atoms to Plotly traces.

        Args:
            config: Configuration for atom visualization appearance.

        Returns:
            List of Plotly Scatter3d traces representing the atoms.
        """
        # Group atoms by element
        element_groups: Dict[str, List[Atom]] = {}
        for atom in self.atoms:
            element_groups.setdefault(atom.element, []).append(atom)

        # Sort elements alphabetically
        sorted_elements = sorted(element_groups.keys())

        traces = []
        for element in sorted_elements:
            atoms = element_groups[element]
            # Convert to numpy arrays for efficient operations
            coords = np.array([[atom.x, atom.y, atom.z] for atom in atoms])
            # Round coordinates to 3 decimal places
            coords = np.round(coords, 3)

            traces.append(
                go.Scatter3d(
                    x=coords[:, 0],
                    y=coords[:, 1],
                    z=coords[:, 2],
                    mode="markers",
                    text=[atom.element for atom in atoms],
                    textposition="top center",
                    marker=dict(
                        size=config["size"],
                        color=self.colors.get(element, "pink"),
                        opacity=config["opacity"],
                        line=dict(color=config["border_color"], width=config["border_width"]),
                    ),
                    hovertext=[
                        f"{atom.index}. {atom.element} ({x:.3f}, {y:.3f}, {z:.3f})"
                        for atom, (x, y, z) in zip(atoms, coords)
                    ],
                    hoverinfo="text",
                    name=element,
                )
            )
        return traces


@dataclass
class BondShape:
    """A class for managing bond shapes in the visualization.

    Handles the conversion of bond data into Plotly visualization traces,
    including grouping bonds by type and applying visual styling.
    """

    bonds: List[Bond]

    def to_trace(self, config: BondConfig = DEFAULT_BOND_CONFIG) -> List[go.Scatter3d]:
        """Convert bonds to Plotly traces.

        Args:
            config: Configuration for bond visualization appearance.

        Returns:
            List of Plotly Scatter3d traces representing the bonds.
        """
        # Group bonds by type
        bond_groups: Dict[str, List[Bond]] = {}
        for bond in self.bonds:
            elements = sorted([bond.atom1.element, bond.atom2.element])
            bond_type = f"{elements[0]}-{elements[1]}"
            bond_groups.setdefault(bond_type, []).append(bond)

        # Sort bond types alphabetically
        sorted_bond_types = sorted(bond_groups.keys())

        traces = []
        for bond_type in sorted_bond_types:
            bonds = bond_groups[bond_type]
            # Create numpy arrays for start and end points
            start_points = np.array([[bond.atom1.x, bond.atom1.y, bond.atom1.z] for bond in bonds])
            end_points = np.array([[bond.atom2.x, bond.atom2.y, bond.atom2.z] for bond in bonds])

            # Round coordinates to 3 decimal places
            start_points = np.round(start_points, 3)
            end_points = np.round(end_points, 3)

            # Calculate midpoints efficiently and round
            midpoints = np.round((start_points + end_points) / 2, 3)

            # Calculate lengths efficiently and round
            lengths = np.round(np.sqrt(np.sum((end_points - start_points) ** 2, axis=1)), 3)
            label_texts = [f"{length:.2f}" for length in lengths]

            # Create bond line coordinates with None separators
            bond_coordinates = np.zeros((len(bonds) * 3, 3))
            bond_coordinates[::3] = start_points
            bond_coordinates[1::3] = end_points
            bond_coordinates[2::3] = np.nan

            # Add bond lines trace (coordinates are already rounded)
            traces.append(
                go.Scatter3d(
                    x=bond_coordinates[:, 0],
                    y=bond_coordinates[:, 1],
                    z=bond_coordinates[:, 2],
                    mode="lines",
                    line=dict(color=config["color"], width=config["width"]),
                    hoverinfo="skip",
                    name=bond_type,
                    legendgroup=bond_type,
                    showlegend=True,
                )
            )

            # Add labels trace (midpoints are already rounded)
            traces.append(
                go.Scatter3d(
                    x=midpoints[:, 0],
                    y=midpoints[:, 1],
                    z=midpoints[:, 2],
                    mode="text",
                    text=label_texts,
                    textposition="top center",
                    hoverinfo="skip",
                    name=f"{bond_type} lengths",
                    legendgroup=bond_type,
                    showlegend=False,
                    textfont=dict(size=10),
                    visible=False,
                )
            )

        return traces


class ButtonFactory:
    """Factory for creating visualization control buttons."""

    def __init__(self, buttons_config: ButtonsConfig) -> None:
        """Initialize the button factory.

        Args:
            buttons_config: Configuration for button appearance and behavior.
        """
        self.config = buttons_config

    def create_atom_buttons(self, atom_traces: List[int]) -> List[Dict[str, Any]]:
        """Create buttons for controlling atom visualization.

        Args:
            atom_traces: List of trace indices corresponding to atom visualizations.

        Returns:
            List of button configurations.
        """
        return [
            {
                "label": "Basic",
                "method": "restyle",
                "args": [{"mode": "markers", "visible": True}, atom_traces],
            },
            {
                "label": "Detailed",
                "method": "restyle",
                "args": [{"mode": "markers+text", "visible": True}, atom_traces],
            },
            {
                "label": "Hide",
                "method": "restyle",
                "args": [{"visible": False}, atom_traces],
            },
        ]

    def create_bond_buttons(self, bond_traces: List[int]) -> List[Dict[str, Any]]:
        """Create buttons for controlling bond visualization.

        Args:
            bond_traces: List of trace indices corresponding to bond visualizations.

        Returns:
            List of button configurations.
        """
        return [
            {
                "label": "Basic",
                "method": "restyle",
                "args": [{"mode": ["lines", "none"], "visible": [True, True]}, bond_traces],
            },
            {
                "label": "Detailed",
                "method": "restyle",
                "args": [{"mode": ["lines", "text"], "visible": [True, True]}, bond_traces],
            },
            {
                "label": "Hide",
                "method": "restyle",
                "args": [{"visible": [False, False]}, bond_traces],
            },
        ]

    def create_grid_buttons(self) -> List[Dict[str, Any]]:
        """Create buttons for controlling grid visualization.

        Returns:
            List of button configurations.
        """

        def make_args(show: bool) -> Dict[str, Any]:
            args: Dict[str, Any] = {}
            for axis in ["xaxis", "yaxis", "zaxis"]:
                args[f"scene.{axis}.showbackground"] = show
                args[f"scene.{axis}.gridcolor"] = "lightgray" if show else "white"
                args[f"scene.{axis}.showgrid"] = show
                args[f"scene.{axis}.backgroundcolor"] = "white" if show else None
            return args

        return [
            {"label": "Show", "method": "relayout", "args": [make_args(True)]},
            {"label": "Hide", "method": "relayout", "args": [make_args(False)]},
        ]

    def create_button_group(
        self, group_name: str, buttons: List[Dict[str, Any]], active_index: int = 0
    ) -> Tuple[Dict[str, Any], Dict[str, Any]]:
        """Create a group of related buttons with associated annotations.

        Args:
            group_name: Name of the button group.
            buttons: List of button configurations.
            active_index: Index of the initially active button.

        Returns:
            Tuple of button group configuration and title annotation.
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
            "y": group_position["label_y"],
            "yref": "paper",
            "align": "left",
            "showarrow": False,
            "font": {"size": appearance["title_font_size"]},  # Updated to use title font size
        }

        return button_group, title_annotation


class VisualizationControls:
    """Manages the interactive controls for the visualization."""

    def __init__(self, n_atoms: int, n_bonds: int, buttons_config: ButtonsConfig) -> None:
        """Initialize visualization controls.

        Args:
            n_atoms: Number of atoms in the visualization.
            n_bonds: Number of bonds in the visualization.
            buttons_config: Configuration for control buttons.
        """
        self.n_atoms = n_atoms
        self.n_bonds = n_bonds
        self.button_factory = ButtonFactory(buttons_config)
        self.atom_traces: List[int] = []
        self.bond_traces: List[int] = []

    def update_trace_indices(self, atom_traces: List[int], bond_traces: List[int]) -> None:
        """Update the indices of atom and bond traces.

        Args:
            atom_traces: List of atom trace indices.
            bond_traces: List of bond trace indices.
        """
        self.atom_traces = atom_traces
        self.bond_traces = bond_traces

    def create_all_controls(self) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
        """Create all visualization controls.

        Returns:
            Tuple of button configurations and annotations.
        """
        buttons = []
        annotations = []

        # Atom controls
        atom_buttons = self.button_factory.create_atom_buttons(self.atom_traces)
        button_group, annotation = self.button_factory.create_button_group("atoms", atom_buttons)
        buttons.append(button_group)
        annotations.append(annotation)

        # Bond controls
        if self.n_bonds > 0:
            bond_buttons = self.button_factory.create_bond_buttons(self.bond_traces)
            button_group, annotation = self.button_factory.create_button_group(
                "bonds", bond_buttons
            )
            buttons.append(button_group)
            annotations.append(annotation)

        # Background controls
        background_buttons = self.button_factory.create_grid_buttons()
        button_group, annotation = self.button_factory.create_button_group(
            "grid", background_buttons
        )
        buttons.append(button_group)
        annotations.append(annotation)

        return buttons, annotations


class VisualizationManager:
    """Manages the overall visualization process."""

    def __init__(self, mol_graph: "MolGraph", config: Optional[VisualizationConfig] = None) -> None:
        """Initialize the visualization manager.

        Args:
            mol_graph: Molecular graph to visualize.
            config: Optional custom visualization configuration.
        """
        self.config_manager = ConfigurationManager(config)
        self.atom_component = AtomShape(atoms=mol_graph.atoms, colors=mol_graph.cpk_colors)
        self.bond_component = BondShape(bonds=mol_graph.bonds)
        self.controls = VisualizationControls(
            n_atoms=len(mol_graph.atoms),
            n_bonds=len(mol_graph.bonds),
            buttons_config=self.config_manager.buttons,
        )

    def _create_axis_config(self) -> Dict[str, Dict[str, Any]]:
        """Create axis configuration for all three axes."""
        axis_config = {
            "showbackground": self.config_manager.scene["showbackground"],
            "showticklabels": self.config_manager.scene["showticklabels"],
            "zeroline": self.config_manager.scene["zeroline"],
            "showspikes": self.config_manager.scene["showspikes"],
            "visible": self.config_manager.scene["axis_visible"],
            "title": {"text": ""},
            "showgrid": self.config_manager.scene["showbackground"],  # Added
            "gridcolor": "lightgray"
            if self.config_manager.scene["showbackground"]
            else "white",  # Added
            "backgroundcolor": "white"
            if self.config_manager.scene["showbackground"]
            else None,  # Modified
        }

        return {axis: axis_config.copy() for axis in ["xaxis", "yaxis", "zaxis"]}

    def create_figure(self) -> go.Figure:
        """Create the Plotly figure for visualization.

        Returns:
            Plotly Figure object containing the molecular visualization.
        """
        fig = go.Figure()

        # Add atom traces
        atom_traces = self.atom_component.to_trace(config=self.config_manager.atoms)
        atom_trace_indices = []
        for trace in atom_traces:
            fig.add_trace(trace)
            atom_trace_indices.append(len(fig.data) - 1)

        # Add bond traces
        bond_traces = self.bond_component.to_trace(config=self.config_manager.bonds)
        bond_trace_indices = []
        for trace in bond_traces:
            fig.add_trace(trace)
            bond_trace_indices.append(len(fig.data) - 1)

        # Update control indices
        self.controls.update_trace_indices(atom_trace_indices, bond_trace_indices)

        # Create buttons and annotations
        buttons, control_annotations = self.controls.create_all_controls()

        # Configure axes
        axis_config = self._create_axis_config()

        # Calculate top margin based on title presence
        top_margin = 40 if self.config_manager.style["title"] else 0

        # Update layout with all configurations
        fig.update_layout(
            title=self.config_manager.style["title"],
            scene={
                **axis_config,
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
