"""Manage and visualize molecular structure traces using Plotly.

This module provides classes for managing and visualizing different types of molecular structure
traces in Plotly visualizations, following a similar pattern to the labels system.
"""

from abc import ABC, abstractmethod
from enum import Enum, auto
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Set, Tuple, Type

import plotly.graph_objects as go

from xyz2graph.layout_config import LAYOUT

from .geometry import Point3D


# Only import types for type checking to avoid circular imports
if TYPE_CHECKING:
    from .graph import MolGraph


class ObjectType(Enum):
    """Enum defining different types of molecular traces."""

    ATOMS = auto()
    BONDS = auto()


class BaseShape(ABC):
    """Abstract base class for structure traces."""

    def __init__(
        self, position: List[Point3D], color: Optional[str] = None, visible: bool = True
    ) -> None:
        """Initialize a base trace."""
        self.position = position
        self.color = color
        self.visible = visible

    def to_plotly_trace(self) -> go.Scatter3d:
        """Convert the trace to a Plotly trace object."""
        return self._create_trace()

    @abstractmethod
    def _create_trace(self) -> go.Scatter3d:
        """Create the specific Plotly trace implementation."""
        pass

    @abstractmethod
    def get_shape_type(self) -> ObjectType:
        """Return the type of this trace."""
        pass


class Atom(BaseShape):
    """Class representing atom traces in molecular structures."""

    @classmethod
    def from_mol_graph(
        cls: Type["Atom"],
        mol_graph: "MolGraph",
        config: Optional[Dict[str, Any]] = None,
    ) -> "Atom":
        """Create an atom trace from molecular graph data."""
        coordinates = [
            Point3D(x=x, y=y, z=z) for x, y, z in zip(mol_graph.x, mol_graph.y, mol_graph.z)
        ]

        return cls(
            position=coordinates,
            elements=mol_graph.elements,
            cpk_colors=mol_graph.cpk_colors,
            default_color=config.get("atom_color", "pink") if config else "pink",
            size=config.get("atom_size", 7) if config else 7,
            opacity=config.get("atom_opacity", 0.8) if config else 0.8,
            border_color=config.get("atom_border_color", "lightgray") if config else "lightgray",
            border_width=config.get("atom_border_width", 2) if config else 2,
        )

    def __init__(
        self,
        position: List[Point3D],
        elements: List[str],
        cpk_colors: Dict[str, str],
        default_color: str = "pink",
        size: int = 7,
        opacity: float = 0.8,
        border_color: str = "lightgray",
        border_width: int = 2,
        visible: bool = True,
    ) -> None:
        """Initialize an atom trace."""
        super().__init__(position, None, visible)
        self.elements = elements
        self.cpk_colors = cpk_colors
        self.default_color = default_color
        self.size = size
        self.opacity = opacity
        self.border_color = border_color
        self.border_width = border_width

    def _create_trace(self) -> go.Scatter3d:
        """Create a Plotly trace for atoms."""
        colors = [self.cpk_colors.get(element, self.default_color) for element in self.elements]

        x = [pos.x for pos in self.position]
        y = [pos.y for pos in self.position]
        z = [pos.z for pos in self.position]

        hover_text = [
            f"{i}. {elem} ({pos.x:.2f}, {pos.y:.2f}, {pos.z:.2f})"
            for i, (elem, pos) in enumerate(zip(self.elements, self.position))
        ]

        return go.Scatter3d(
            x=x,
            y=y,
            z=z,
            mode="markers",
            marker=dict(
                color=colors,
                size=self.size,
                opacity=self.opacity,
                line=dict(
                    color=self.border_color,
                    width=self.border_width,
                ),
                symbol="circle",
            ),
            text=self.elements,
            hovertext=hover_text,
            hoverinfo="text",
            name="atoms",
            visible=self.visible,
        )

    def get_shape_type(self) -> ObjectType:
        """Return the type of this trace."""
        return ObjectType.ATOMS


class Bond(BaseShape):
    """Class representing bond traces in molecular structures."""

    @classmethod
    def from_mol_graph(
        cls: Type["Bond"],
        mol_graph: "MolGraph",
        config: Optional[Dict[str, Any]] = None,
    ) -> "Bond":
        """Create a bond trace from molecular graph data."""
        coordinates = [
            Point3D(x=x, y=y, z=z) for x, y, z in zip(mol_graph.x, mol_graph.y, mol_graph.z)
        ]

        return cls(
            position=coordinates,
            adjacency_list=mol_graph.adj_list,
            color=config.get("bond_color", "grey") if config else "grey",
            width=config.get("bond_width", 2) if config else 2,
        )

    def __init__(
        self,
        position: List[Point3D],
        adjacency_list: Dict[int, Set[int]],
        color: str = "grey",
        width: int = 2,
        visible: bool = True,
    ) -> None:
        """Initialize a bond trace."""
        super().__init__(position, color, visible)
        self.adjacency_list = adjacency_list
        self.width = width

    def _create_trace(self) -> go.Scatter3d:
        """Create a Plotly trace for bonds."""
        xs, ys, zs = [], [], []
        processed_bonds = set()

        for atom1, neighbors in self.adjacency_list.items():
            for atom2 in neighbors:
                bond = frozenset([atom1, atom2])
                if bond in processed_bonds:
                    continue

                processed_bonds.add(bond)

                coord1 = self.position[atom1]
                coord2 = self.position[atom2]

                xs.extend([coord1.x, coord2.x, None])
                ys.extend([coord1.y, coord2.y, None])
                zs.extend([coord1.z, coord2.z, None])

        return go.Scatter3d(
            x=xs,
            y=ys,
            z=zs,
            mode="lines",
            line=dict(
                color=self.color,
                width=self.width,
            ),
            hoverinfo="none",
            name="bonds",
            visible=self.visible,
        )

    def get_shape_type(self) -> ObjectType:
        """Return the type of this trace."""
        return ObjectType.BONDS


class ButtonFactory:
    """Factory class for creating Plotly button configurations."""

    @staticmethod
    def create_toggle_button(label: str, trace_index: int, visible: bool) -> Dict[str, Any]:
        """Create a single toggle button configuration.

        Args:
            label: Text to display on the button
            trace_indices: Indices of traces to affect
            visible: Visibility state to set

        Returns:
            Dictionary containing Plotly button configuration
        """
        return {
            "label": label,
            "method": "restyle",
            "args": [{"visible": visible}, [trace_index]],
        }


class TraceManager:
    """Manager class for molecular structure traces."""

    def __init__(self) -> None:
        """Initialize an empty trace manager."""
        self._traces: List[BaseShape] = []

    def add_trace(self, trace: BaseShape) -> None:
        """Add a new trace to the manager.

        Args:
            trace: BaseTrace instance to add
        """
        self._traces.append(trace)

    def clear(self) -> None:
        """Remove all traces from the manager."""
        self._traces.clear()

    def get_traces_by_type(self, trace_type: ObjectType) -> List[BaseShape]:
        """Get all traces of a specific type.

        Args:
            trace_type: Type of traces to retrieve

        Returns:
            List of traces matching the specified type
        """
        return [trace for trace in self._traces if trace.get_shape_type() == trace_type]

    def to_plotly_traces(self) -> List[go.Scatter3d]:
        """Convert all traces to Plotly format.

        Returns:
            List of Plotly trace objects
        """
        return [trace.to_plotly_trace() for trace in self._traces]

    def _get_trace_indices_by_type(self, trace_type: ObjectType) -> List[int]:
        """Get indices of all traces of a specific type.

        Args:
            trace_type: Type of traces to find

        Returns:
            List of indices where traces match the specified type
        """
        return [i for i, trace in enumerate(self._traces) if trace.get_shape_type() == trace_type]

    def create_buttons(self) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
        """Create Plotly buttons for toggling trace visibility.

        Creates separate show/hide buttons for atoms and bonds.

        Returns:
            Tuple containing:
                - List of button group configurations
                - List of group title annotations
        """

        def create_button_group(
            y_position: float, buttons_config: List[Dict[str, Any]], active_index: int = 0
        ) -> Dict[str, Any]:
            """Helper function to create a button group configuration."""
            return {
                "buttons": buttons_config,
                "direction": "right",
                "pad": {"r": 10, "t": 10},
                "showactive": True,
                "active": active_index,
                "x": LAYOUT["common"]["button_x"],
                "y": y_position,
                "xanchor": "left",
                "yanchor": "top",
                "type": "buttons",
            }

        def create_group_annotation(title: str, y_position: float) -> Dict[str, Any]:
            """Helper function to create a group annotation."""
            return {
                "text": title,
                "x": LAYOUT["common"]["text_x"],
                "xref": "paper",
                "y": y_position,
                "yref": "paper",
                "align": "left",
                "showarrow": False,
            }

        buttons = []
        annotations = []

        # Get indices for atoms and bonds
        atom_indices = self._get_trace_indices_by_type(ObjectType.ATOMS)
        bond_indices = self._get_trace_indices_by_type(ObjectType.BONDS)

        if atom_indices:
            atom_buttons = [
                ButtonFactory.create_toggle_button("Show", atom_indices[0], True),
                ButtonFactory.create_toggle_button("Hide", atom_indices[0], False),
            ]
            buttons.append(
                create_button_group(
                    LAYOUT["atom_traces"]["button_y"],
                    atom_buttons,
                )
            )

            annotations.append(
                create_group_annotation(
                    LAYOUT["atom_traces"]["title"],
                    LAYOUT["atom_traces"]["text_y"],
                )
            )

        if bond_indices:
            bond_buttons = [
                ButtonFactory.create_toggle_button("Show", bond_indices[0], True),
                ButtonFactory.create_toggle_button("Hide", bond_indices[0], False),
            ]
            buttons.append(
                create_button_group(
                    LAYOUT["bond_traces"]["button_y"],
                    bond_buttons,
                )
            )
            annotations.append(
                create_group_annotation(
                    LAYOUT["bond_traces"]["title"],
                    LAYOUT["bond_traces"]["text_y"],
                )
            )

        return buttons, annotations

    def update_figure(self, fig: go.Figure) -> None:
        """Add traces to figure.

        Args:
            fig: Plotly figure to update
        """
        for trace in self.to_plotly_traces():
            fig.add_trace(trace)

    def get_menu_items(self) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
        """Get button groups and annotations for trace controls.

        Returns:
            Tuple containing:
                - List of button group configurations
                - List of group title annotations
        """
        return self.create_buttons()
