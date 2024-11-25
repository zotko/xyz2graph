"""Visualization module for molecular structures using Plotly.

This module provides classes and utilities for creating interactive 3D
visualizations of molecular structures. It includes support for atoms,
bonds, labels, and interactive controls.
"""

from dataclasses import dataclass
from textwrap import shorten
from typing import (
    TYPE_CHECKING,
    Any,
    Dict,
    List,
    Optional,
    Tuple,
    TypedDict,
)

import numpy as np
import plotly.graph_objects as go

from .molecule import Atom, Bond


if TYPE_CHECKING:
    from .graph import MolGraph


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
    """Position configuration for a button group."""

    button_x: float
    button_y: float
    label_x: float
    label_y: float
    title: str


class ButtonsConfig(TypedDict):
    """Complete button configuration."""

    appearance: ButtonAppearance
    groups: Dict[str, ButtonPosition]


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
    text_size: int


class LegendConfig(TypedDict):
    """Legend display configuration."""

    x: float
    y: float
    xanchor: str
    yanchor: str
    itemsizing: str


class WatermarkConfig(TypedDict):
    """Watermark configuration."""

    text: str
    x: float
    y: float
    font_size: int
    font_color: str
    show: bool


class FormulaConfig(TypedDict):
    """Formula annotation configuration."""

    x: float
    y: float
    font_size: int
    font_color: str
    xanchor: str
    yanchor: str


class CommentConfig(TypedDict):
    """Comment annotation configuration."""

    x: float
    y: float
    font_size: int
    font_color: str
    xanchor: str
    yanchor: str
    max_length: int


class InformationConfig(TypedDict):
    """Information display configuration."""

    formula: FormulaConfig
    comment: CommentConfig


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
    watermark: WatermarkConfig
    information: InformationConfig
    coordinate_round_digits: int
    markdown_round_digits: int


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

DEFAULT_BUTTON_X = 0.04
DEFAULT_LABEL_X = 0.0

DEFAULT_BUTTON_GROUPS: Dict[str, ButtonPosition] = {
    "atoms": {
        "button_x": DEFAULT_BUTTON_X,
        "button_y": 1.0,
        "label_x": DEFAULT_LABEL_X,
        "label_y": 0.99,
        "title": "Atoms",
    },
    "bonds": {
        "button_x": DEFAULT_BUTTON_X,
        "button_y": 0.95,
        "label_x": DEFAULT_LABEL_X,
        "label_y": 0.94,
        "title": "Bonds",
    },
    "grid": {
        "button_x": DEFAULT_BUTTON_X,
        "button_y": 0.90,
        "label_x": DEFAULT_LABEL_X,
        "label_y": 0.89,
        "title": "Grid",
    },
}

DEFAULT_BUTTONS = ButtonsConfig(
    appearance=DEFAULT_BUTTON_APPEARANCE,
    groups=DEFAULT_BUTTON_GROUPS,
)

DEFAULT_ATOM_CONFIG = AtomConfig(
    size=6,
    opacity=0.9,
    border_color="black",
    border_width=4,
    hover_font_size=12,
)

DEFAULT_BOND_CONFIG = BondConfig(
    color="grey",
    width=2,
    text_size=10,
)

DEFAULT_LEGEND_CONFIG = LegendConfig(
    x=0.99, y=0.96, xanchor="right", yanchor="top", itemsizing="constant"
)

DEFAULT_WATERMARK_CONFIG = WatermarkConfig(
    text=(
        'created with <a href="https://zotko.github.io/xyz2graph/"'
        'style="color:#404040;">xyz2graph</a>'
    ),
    x=0.99,
    y=0.01,
    font_size=10,
    font_color="gray",
    show=True,
)

DEFAULT_FORMULA_CONFIG = FormulaConfig(
    x=0.00,
    y=0.03,
    font_size=12,
    font_color="black",
    xanchor="left",
    yanchor="bottom",
)

DEFAULT_COMMENT_CONFIG = CommentConfig(
    x=0.00,
    y=0.01,
    font_size=10,
    font_color="black",
    xanchor="left",
    yanchor="bottom",
    max_length=100,
)

DEFAULT_INFORMATION_CONFIG = InformationConfig(
    formula=DEFAULT_FORMULA_CONFIG,
    comment=DEFAULT_COMMENT_CONFIG,
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

DEFAULT_STYLE = StyleConfig(
    atoms=DEFAULT_ATOM_CONFIG,
    bonds=DEFAULT_BOND_CONFIG,
    watermark=DEFAULT_WATERMARK_CONFIG,
    information=DEFAULT_INFORMATION_CONFIG,
    coordinate_round_digits=3,
    markdown_round_digits=2,
)

DEFAULT_CONFIG = VisualizationConfig(
    style=DEFAULT_STYLE, scene=DEFAULT_SCENE_CONFIG, buttons=DEFAULT_BUTTONS
)


class ConfigurationManager:
    """Manages and provides access to visualization configuration components.

    This class serves as a central configuration manager that provides structured
    access to different aspects of the visualization configuration including:
    - Style settings (atoms, bonds, information, watermark)
    - Scene settings (background, axes, legend)
    - Button controls configuration
    - Atom and bond-specific visualization parameters

    The configuration is initialized with either custom settings or defaults,
    and provides property-based access to each configuration component.
    """

    def __init__(self, config: Optional[VisualizationConfig] = None) -> None:
        """Initialize configuration manager.

        Args:
            config: Custom visualization configuration. Uses DEFAULT_CONFIG if None.
        """
        self.config = config or DEFAULT_CONFIG

    @property
    def style(self) -> StyleConfig:
        """Get the style configuration settings.

        Returns:
            StyleConfig: Configuration for visual styling including atoms,
                bonds, title, and watermark settings.
        """
        return self.config["style"]

    @property
    def scene(self) -> SceneConfig:
        """Get the scene configuration settings.

        Returns:
            SceneConfig: Configuration for 3D scene including background,
                axes, and legend settings.
        """
        return self.config["scene"]

    @property
    def buttons(self) -> ButtonsConfig:
        """Get the buttons configuration settings.

        Returns:
            ButtonsConfig: Configuration for interactive control buttons
                including appearance and positioning.
        """
        return self.config["buttons"]

    @property
    def atoms(self) -> AtomConfig:
        """Get the atom visualization configuration.

        Returns:
            AtomConfig: Configuration specific to atom visualization
                including size, opacity, and border settings.
        """
        return self.style["atoms"]

    @property
    def bonds(self) -> BondConfig:
        """Get the bond visualization configuration.

        Returns:
            BondConfig: Configuration specific to bond visualization
                including color, width, and text settings.
        """
        return self.style["bonds"]

    @property
    def watermark(self) -> WatermarkConfig:
        """Get the watermark configuration.

        Returns:
            WatermarkConfig: Configuration for the visualization watermark
                including text, position, and appearance settings.
        """
        return self.style["watermark"]


@dataclass
class AtomShape:
    """A class for managing atom shapes in the visualization."""

    atoms: List[Atom]
    colors: Dict[str, str]

    def to_trace(self, config: AtomConfig, coordinate_round_digits: int) -> List[go.Scatter3d]:
        """Convert atoms to Plotly visualization traces.

        This method processes atoms into Plotly Scatter3d traces for visualization,
        grouping atoms by element and applying the specified visual styling.
        Each element gets its own trace with appropriate coloring and hover text.

        Args:
            config (AtomConfig): Configuration for atom appearance including size,
                opacity, and border settings.
            coordinate_round_digits (int): Number of decimal places to round
                coordinate values.

        Returns:
            List[go.Scatter3d]: List of Plotly traces, one for each unique
                element present in the atoms list. Each trace contains the 3D
                coordinates, styling, and hover information for all atoms of
                that element.
        """
        element_groups: Dict[str, List[Atom]] = {}
        for atom in self.atoms:
            element_groups.setdefault(atom.element, []).append(atom)

        sorted_elements = sorted(element_groups.keys())

        traces = []
        for element in sorted_elements:
            atoms = element_groups[element]
            coords = np.array([[atom.x, atom.y, atom.z] for atom in atoms])
            coords = np.round(coords, coordinate_round_digits)

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
                        f"{atom.index}. {atom.element} ({x:.{coordinate_round_digits}f}, "
                        f"{y:.{coordinate_round_digits}f}, {z:.{coordinate_round_digits}f})"
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

    def to_trace(
        self, config: BondConfig, coordinate_round_digits: int, markdown_round_digits: int
    ) -> List[go.Scatter3d]:
        """Convert bonds to Plotly visualization traces.

        Processes molecular bonds into Plotly Scatter3d traces, grouping bonds by
        type (e.g., C-C, C-O) and creating both line traces for the bonds and
        text traces for the bond lengths.

        Args:
            config (BondConfig): Configuration for bond visualization appearance,
                including color, width, and text settings.
            coordinate_round_digits (int): Number of decimal places to round
                coordinate values in the visualization.
            markdown_round_digits (int): Number of decimal places to round
                bond length values displayed in the visualization.

        Returns:
            List[go.Scatter3d]: List of Plotly traces where:
                - Even-indexed traces represent bond lines
                - Odd-indexed traces represent bond length labels
                Each bond type (e.g., C-C, C-O) gets its own pair of traces.
        """
        # Group bonds by type
        bond_groups: Dict[str, List[Bond]] = {}
        for bond in self.bonds:
            elements = sorted([bond.atom1.element, bond.atom2.element])
            bond_type = f"{elements[0]}-{elements[1]}"
            bond_groups.setdefault(bond_type, []).append(bond)

        sorted_bond_types = sorted(bond_groups.keys())

        traces = []
        for bond_type in sorted_bond_types:
            bonds = bond_groups[bond_type]

            start_points = np.array([[bond.atom1.x, bond.atom1.y, bond.atom1.z] for bond in bonds])
            end_points = np.array([[bond.atom2.x, bond.atom2.y, bond.atom2.z] for bond in bonds])
            midpoints = np.round((start_points + end_points) / 2, coordinate_round_digits)

            lengths = np.round(
                np.sqrt(np.sum((end_points - start_points) ** 2, axis=1)), markdown_round_digits
            )

            bond_coordinates = np.zeros((len(bonds) * 3, 3))
            bond_coordinates[::3] = np.round(start_points, coordinate_round_digits)
            bond_coordinates[1::3] = np.round(end_points, coordinate_round_digits)
            bond_coordinates[2::3] = np.nan

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

            traces.append(
                go.Scatter3d(
                    x=midpoints[:, 0],
                    y=midpoints[:, 1],
                    z=midpoints[:, 2],
                    mode="text",
                    text=[f"{length:.{markdown_round_digits}f}" for length in lengths],
                    textposition="top center",
                    hoverinfo="skip",
                    name=f"{bond_type} lengths",
                    legendgroup=bond_type,
                    showlegend=False,
                    textfont=dict(size=config["text_size"]),
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
        group_config = self.config["groups"][group_name]
        appearance = self.config["appearance"]

        button_group = {
            "buttons": buttons,
            "direction": appearance["direction"],
            "pad": {"r": appearance["pad_r"], "t": appearance["pad_t"]},
            "showactive": appearance["showactive"],
            "active": active_index,
            "x": group_config["button_x"],
            "y": group_config["button_y"],
            "xanchor": appearance["xanchor"],
            "yanchor": appearance["yanchor"],
            "type": "buttons",
            "font": {"size": appearance["font_size"]},
        }

        title_annotation = {
            "text": group_config["title"],
            "x": group_config["label_x"],
            "y": group_config["label_y"],
            "yref": "paper",
            "align": "left",
            "showarrow": False,
            "font": {"size": appearance["title_font_size"]},
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
        self.mol_graph = mol_graph
        self.config_manager = ConfigurationManager(config)
        self.atom_component = AtomShape(atoms=mol_graph.atoms, colors=mol_graph.cpk_colors)
        self.bond_component = BondShape(bonds=mol_graph.bonds)
        self.controls = VisualizationControls(
            n_atoms=len(mol_graph.atoms),
            n_bonds=len(mol_graph.bonds),
            buttons_config=self.config_manager.buttons,
        )

    def _create_axis_config(self) -> Dict[str, Dict[str, Any]]:
        """Create the configuration dictionary for scene axes.

        Generates a configuration dictionary for the x, y, and z axes of the 3D
        scene, applying settings from the scene configuration including grid,
        background, and visibility settings.

        Returns:
            Dict[str, Dict[str, Any]]: Dictionary containing identical
                configuration for all three axes (xaxis, yaxis, zaxis) with
                settings for background, grid, ticks, and visibility.
        """
        axis_config = {
            "showbackground": self.config_manager.scene["showbackground"],
            "showticklabels": self.config_manager.scene["showticklabels"],
            "zeroline": self.config_manager.scene["zeroline"],
            "showspikes": self.config_manager.scene["showspikes"],
            "visible": self.config_manager.scene["axis_visible"],
            "title": {"text": ""},
            "showgrid": self.config_manager.scene["showbackground"],
            "gridcolor": "lightgray"
            if self.config_manager.scene["showbackground"]
            else "white",  # Added
            "backgroundcolor": "white"
            if self.config_manager.scene["showbackground"]
            else None,  # Modified
        }

        return {axis: axis_config.copy() for axis in ["xaxis", "yaxis", "zaxis"]}

    def _create_annotation(
        self,
        text: str,
        x: float,
        y: float,
        font_size: int,
        font_color: str,
        xanchor: str = "left",
        yanchor: str = "bottom",
    ) -> Dict[str, Any]:
        """Create a standardized annotation dictionary.

        Args:
            text: Text content for the annotation
            x: X-position (0-1 range)
            y: Y-position (0-1 range)
            font_size: Size of the font
            font_color: Color of the font
            xanchor: Horizontal anchor point ("left" or "right")
            yanchor: Vertical anchor point ("bottom" or "top")

        Returns:
            Dict[str, Any]: Annotation configuration dictionary
        """
        return {
            "text": text,
            "xref": "paper",
            "yref": "paper",
            "x": x,
            "y": y,
            "showarrow": False,
            "font": {
                "size": font_size,
                "color": font_color,
            },
            "xanchor": xanchor,
            "yanchor": yanchor,
        }

    def create_figure(self) -> go.Figure:
        """Create the Plotly figure for visualization.

        Returns:
            Plotly Figure object containing the molecular visualization.
        """
        fig = go.Figure()

        # Add atom traces
        atom_traces = self.atom_component.to_trace(
            config=self.config_manager.atoms,
            coordinate_round_digits=self.config_manager.style["coordinate_round_digits"],
        )
        atom_trace_indices = []
        for trace in atom_traces:
            fig.add_trace(trace)
            atom_trace_indices.append(len(fig.data) - 1)

        # Add bond traces
        bond_traces = self.bond_component.to_trace(
            config=self.config_manager.bonds,
            coordinate_round_digits=self.config_manager.style["coordinate_round_digits"],
            markdown_round_digits=self.config_manager.style["markdown_round_digits"],
        )
        bond_trace_indices = []
        for trace in bond_traces:
            fig.add_trace(trace)
            bond_trace_indices.append(len(fig.data) - 1)

        # Update control indices
        self.controls.update_trace_indices(atom_trace_indices, bond_trace_indices)

        # Create buttons and control annotations
        buttons, control_annotations = self.controls.create_all_controls()

        # Create information annotations
        information_annotations = []
        info_config = self.config_manager.style["information"]

        # Add formula annotation if present
        formula = self.mol_graph.formula()
        if formula:
            formula_config = info_config["formula"]
            information_annotations.append(
                self._create_annotation(
                    text=formula,
                    x=formula_config["x"],
                    y=formula_config["y"],
                    font_size=formula_config["font_size"],
                    font_color=formula_config["font_color"],
                    xanchor=formula_config["xanchor"],
                    yanchor=formula_config["yanchor"],
                )
            )

        # Add comment annotation if present
        comment = shorten(
            self.mol_graph.comment.strip(),
            width=info_config["comment"]["max_length"],
            placeholder="...",
        )
        if comment:
            comment_config = info_config["comment"]
            information_annotations.append(
                self._create_annotation(
                    text=comment,
                    x=comment_config["x"],
                    y=comment_config["y"],
                    font_size=comment_config["font_size"],
                    font_color=comment_config["font_color"],
                    xanchor=comment_config["xanchor"],
                    yanchor=comment_config["yanchor"],
                )
            )

        # Add watermark if enabled
        watermark_config = self.config_manager.watermark
        if watermark_config["show"]:
            information_annotations.append(
                self._create_annotation(
                    text=watermark_config["text"],
                    x=watermark_config["x"],
                    y=watermark_config["y"],
                    font_size=watermark_config["font_size"],
                    font_color=watermark_config["font_color"],
                    xanchor="right",
                    yanchor="bottom",
                )
            )

        # Combine all annotations
        all_annotations = control_annotations + information_annotations

        # Update layout with all configurations
        fig.update_layout(
            scene={
                **self._create_axis_config(),
            },
            updatemenus=buttons,
            annotations=all_annotations,
            margin=dict(r=0, l=0, b=0, t=0),
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
