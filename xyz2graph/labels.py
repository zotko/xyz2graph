"""Manage and visualize molecular structure labels using Plotly.

This module provides classes and functions for managing and visualizing molecular structure
labels in Plotly visualizations.

Classes:
    LabelType: Enum defining different types of molecular labels.
    BaseLabel: Abstract base class for structure labels.
    AtomLabel: Class representing atom labels in molecular structures.
    BondLabel: Class representing bond labels in molecular structures.
    ButtonFactory: Factory class for creating Plotly button configurations.
    MolecularLabelManager: Manager class for molecular structure labels.
"""

from abc import ABC, abstractmethod
from enum import Enum, auto
from typing import Any, Dict, List, Optional, Tuple, Type, TypedDict

import plotly.graph_objects as go

from .geometry import Point3D


# Default layout configuration for buttons and annotations
class _CommonLayout(TypedDict):
    button_x: float
    text_x: float


class _LabelGroupLayout(TypedDict):
    button_y: float
    text_y: float
    title: str


class _LayoutConfig(TypedDict):
    atom_labels: _LabelGroupLayout
    bond_lengths: _LabelGroupLayout
    common: _CommonLayout


LAYOUT: _LayoutConfig = {
    "atom_labels": {"button_y": 1.12, "text_y": 1.1, "title": "Atom Labels"},
    "bond_lengths": {"button_y": 1.065, "text_y": 1.045, "title": "Bond Lengths"},
    "common": {
        "button_x": 0.07,
        "text_x": 0,
    },
}


class LabelType(Enum):
    """Enum defining different types of molecular labels.

    Members:
        ATOM_ELEMENT: Label showing atomic element symbol
        ATOM_INDEX: Label showing atom index
        BOND_LENGTH: Label showing bond length
    """

    ATOM_ELEMENT = auto()
    ATOM_INDEX = auto()
    BOND_LENGTH = auto()


class BaseLabel(ABC):
    """Abstract base class for structure labels.

    Attributes:
        position: 3D coordinates where the label should be placed
        text: Content to be displayed in the label
        color: Optional color for the label text
        offset: Vertical offset from the label position in pixels
    """

    def __init__(
        self, position: Point3D, text: str, color: Optional[str] = None, offset: int = 15
    ) -> None:
        """Initialize a base label.

        Args:
            position: 3D coordinates where the label should be placed
            text: Content to be displayed in the label
            color: Optional color for the label text. Defaults to None.
            offset: Vertical offset from the label position in pixels. Defaults to 15.
        """
        self.position = position
        self.text = text
        self.color = color
        self.offset = offset

    def to_plotly_annotation(self, visible: bool = False) -> Dict[str, Any]:
        """Convert the label to a Plotly annotation dictionary.

        Args:
            visible: Whether the annotation should be visible. Defaults to False.

        Returns:
            Dictionary containing Plotly annotation parameters
        """
        return {
            "text": self.text,
            "x": self.position.x,
            "y": self.position.y,
            "z": self.position.z,
            "showarrow": False,
            "yshift": self.offset,
            "visible": visible,
            **({"font": {"color": self.color}} if self.color else {}),
        }

    @abstractmethod
    def get_label_type(self) -> LabelType:
        """Return the type of this label.

        Returns:
            LabelType indicating the category of this label
        """
        pass


class AtomLabel(BaseLabel):
    """Class representing atom labels in molecular structures."""

    def __init__(
        self,
        position: Point3D,
        text: str,
        atom_index: int,
        element: str,
        label_type: LabelType,
        color: Optional[str] = None,
        offset: int = 15,
    ) -> None:
        """Initialize an atom label.

        Args:
            position: 3D coordinates where the label should be placed
            text: Content to be displayed in the label
            atom_index: Index of the atom being labeled
            element: Chemical element symbol of the atom
            label_type: Type of atom label (element symbol or index)
            color: Optional color for the label text. Defaults to None.
            offset: Vertical offset from the label position in pixels. Defaults to 15.
        """
        super().__init__(position, text, color, offset)
        self.atom_index = atom_index
        self.element = element
        self._label_type = label_type

    def get_label_type(self) -> LabelType:
        """Return the type of this atom label.

        Returns:
            LabelType.ATOM_ELEMENT or LabelType.ATOM_INDEX
        """
        return self._label_type

    @classmethod
    def from_atom(
        cls: Type["AtomLabel"],
        atom_index: int,
        element: str,
        position: Point3D,
        label_type: LabelType,
        color: Optional[str] = None,
        offset: int = 15,
    ) -> "AtomLabel":
        """Create an atom label from atom data.

        Args:
            atom_index: Index of the atom being labeled
            element: Chemical element symbol of the atom
            position: 3D coordinates where the label should be placed
            label_type: Type of atom label (element symbol or index)
            color: Optional color for the label text. Defaults to None.
            offset: Vertical offset from the label position in pixels. Defaults to 15.

        Returns:
            New AtomLabel instance
        """
        text = str(atom_index) if label_type == LabelType.ATOM_INDEX else element
        return cls(
            position=position,
            text=text,
            atom_index=atom_index,
            element=element,
            label_type=label_type,
            color=color,
            offset=offset,
        )


class BondLabel(BaseLabel):
    """Class representing bond labels in molecular structures."""

    def __init__(
        self,
        position: Point3D,
        bond_length: float,
        atom1_index: int,
        atom2_index: int,
        color: Optional[str] = None,
        offset: int = 15,
    ) -> None:
        """Initialize a bond label.

        Args:
            position: 3D coordinates where the label should be placed
            bond_length: Length of the bond in Angstroms
            atom1_index: Index of the first atom in the bond
            atom2_index: Index of the second atom in the bond
            color: Optional color for the label text. Defaults to None.
            offset: Vertical offset from the label position in pixels. Defaults to 15.
        """
        super().__init__(position, f"{bond_length:.2f}", color, offset)
        self.atom1_index = atom1_index
        self.atom2_index = atom2_index
        self.bond_length = bond_length

    def get_label_type(self) -> LabelType:
        """Return the type of this bond label.

        Returns:
            LabelType.BOND_LENGTH
        """
        return LabelType.BOND_LENGTH

    @classmethod
    def from_bond(
        cls: Type["BondLabel"],
        atom1_index: int,
        atom2_index: int,
        atom1_pos: Point3D,
        atom2_pos: Point3D,
        bond_length: float,
        color: Optional[str] = None,
        offset: int = 15,
    ) -> "BondLabel":
        """Create a bond label from bond data.

        Args:
            atom1_index: Index of the first atom in the bond
            atom2_index: Index of the second atom in the bond
            atom1_pos: 3D coordinates of the first atom
            atom2_pos: 3D coordinates of the second atom
            bond_length: Length of the bond in Angstroms
            color: Optional color for the label text. Defaults to None.
            offset: Vertical offset from the label position in pixels. Defaults to 15.

        Returns:
            New BondLabel instance
        """
        position = Point3D.midpoint(atom1_pos, atom2_pos)
        return cls(
            position=position,
            bond_length=bond_length,
            atom1_index=atom1_index,
            atom2_index=atom2_index,
            color=color,
            offset=offset,
        )


class ButtonFactory:
    """Factory class for creating Plotly button configurations."""

    @staticmethod
    def create_toggle_button(
        label: str, visible_indices: List[int], all_indices: List[int]
    ) -> Dict[str, Any]:
        """Create a single toggle button configuration.

        Args:
            label: Text to display on the button
            visible_indices: Indices of annotations to make visible
            all_indices: All possible annotation indices in the figure

        Returns:
            Dictionary containing Plotly button configuration
        """
        return {
            "label": label,
            "method": "update",
            "args": [
                {},
                {
                    **{
                        f"scene.annotations[{i}].visible": (i in visible_indices)
                        for i in all_indices
                    }
                },
            ],
        }


class MolecularLabelManager:
    """Manager class for molecular structure labels.

    Handles creation, organization, and visualization of molecular structure labels in Plotly
    figures. Supports atom labels (element symbols and indices) and bond length labels with
    toggle buttons for controlling visibility.
    """

    def __init__(self) -> None:
        """Initialize an empty label manager."""
        self._labels: List[BaseLabel] = []

    def add_label(self, label: BaseLabel) -> None:
        """Add a new label to the manager.

        Args:
            label: BaseLabel instance to add
        """
        self._labels.append(label)

    def clear(self) -> None:
        """Remove all labels from the manager."""
        self._labels.clear()

    def get_labels_by_type(self, label_type: LabelType) -> List[BaseLabel]:
        """Get all labels of a specific type.

        Args:
            label_type: Type of labels to retrieve

        Returns:
            List of labels matching the specified type
        """
        return [label for label in self._labels if label.get_label_type() == label_type]

    def to_plotly_annotations(self, visible: bool = False) -> List[Dict[str, Any]]:
        """Convert all labels to Plotly annotation format.

        Args:
            visible: Whether the annotations should be visible. Defaults to False.

        Returns:
            List of Plotly annotation dictionaries
        """
        return [label.to_plotly_annotation(visible=visible) for label in self._labels]

    def _get_label_indices_by_type(self, label_type: LabelType) -> List[int]:
        """Get indices of all labels of a specific type.

        Args:
            label_type: Type of labels to find

        Returns:
            List of indices where labels match the specified type
        """
        return [i for i, label in enumerate(self._labels) if label.get_label_type() == label_type]

    def create_buttons(self) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
        """Create Plotly buttons for toggling label visibility.

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

        # Handle Atom Labels
        element_indices = self._get_label_indices_by_type(LabelType.ATOM_ELEMENT)
        index_indices = self._get_label_indices_by_type(LabelType.ATOM_INDEX)

        if element_indices and index_indices:
            all_atom_indices = element_indices + index_indices
            atom_buttons = [
                ButtonFactory.create_toggle_button("Element", element_indices, all_atom_indices),
                ButtonFactory.create_toggle_button("ID", index_indices, all_atom_indices),
                ButtonFactory.create_toggle_button("Hide", [], all_atom_indices),
            ]

            buttons.append(
                create_button_group(LAYOUT["atom_labels"]["button_y"], atom_buttons, active_index=2)
            )

            annotations.append(
                create_group_annotation(
                    LAYOUT["atom_labels"]["title"], LAYOUT["atom_labels"]["text_y"]
                )
            )

        # Handle Bond Lengths
        bond_indices = self._get_label_indices_by_type(LabelType.BOND_LENGTH)

        if bond_indices:
            bond_buttons = [
                ButtonFactory.create_toggle_button("Show", bond_indices, bond_indices),
                ButtonFactory.create_toggle_button("Hide", [], bond_indices),
            ]

            buttons.append(
                create_button_group(
                    LAYOUT["bond_lengths"]["button_y"],
                    bond_buttons,
                    active_index=1,
                )
            )

            annotations.append(
                create_group_annotation(
                    LAYOUT["bond_lengths"]["title"], LAYOUT["bond_lengths"]["text_y"]
                )
            )

        return buttons, annotations

    def update_figure_menu(self, fig: go.Figure) -> None:
        """Update figure layout with toggle buttons and initial annotations.

        Sets initial visibility of labels (atom elements visible by default) and adds
        button groups for controlling label visibility.

        Args:
            fig: Plotly figure to update
        """
        annotations = self.to_plotly_annotations(visible=False)

        fig.update_layout(scene=dict(annotations=annotations))

        # Add button groups and their annotations
        button_groups, group_annotations = self.create_buttons()
        if button_groups:
            fig.update_layout(updatemenus=button_groups, annotations=group_annotations)
