"""Manage and visualize molecular structure labels using Plotly.

This module defines classes and functions for managing and visualizing molecular
structure labels using Plotly.

Classes:
    LabelType(Enum): Enum defining different types of molecular labels.

    BaseLabel(ABC): Abstract base class for structure labels.

Methods:
            to_plotly_annotation: Convert label to Plotly annotation format.
            get_label_type: Return the type of label.

    AtomLabel(BaseLabel): Class representing atom labels in molecular structure.

Methods:
            from_atom: Create an atom label from atom data.
            to_plotly_annotation: Convert label to Plotly annotation format.
            get_label_type: Return the type of label.

    BondLabel(BaseLabel): Class representing bond labels in molecular structure.

Methods:
            from_bond: Create a bond label from bond data.
            to_plotly_annotation: Convert label to Plotly annotation format.
            get_label_type: Return the type of label.

    PlotlyButtonFactory: Factory class for creating Plotly button configurations.

Methods:
            create_toggle_button: Create a single toggle button configuration.

    MolecularLabelManager: Manager class for molecular structure labels.

Methods:
            add_label: Add a new label.
            clear: Remove all labels.
            get_labels_by_type: Get labels of a specific type.
            to_plotly_annotations: Convert all labels to Plotly annotation format.
            create_toggle_buttons: Create Plotly buttons for toggling label visibility in groups.
            update_figure_menu: Update figure layout with grouped toggle buttons
                and initial annotations.

        Private Methods:
            _get_label_indices_by_type: Get indices of labels of a specific type.
"""

from abc import ABC, abstractmethod
from enum import Enum, auto
from typing import Any, Dict, List, Optional, Tuple, Type

import plotly.graph_objects as go

from .geometry import Point3D


class LabelType(Enum):
    """Enum defining different types of molecular labels."""

    ATOM_ELEMENT = auto()
    ATOM_INDEX = auto()
    BOND_LENGTH = auto()


class BaseLabel(ABC):
    """Abstract base class for structure labels.

    Attributes:
        position (Point3D): The 3D position of the label.
        text (str): The text content of the label.
        color (Optional[str]): The color of the label text. Defaults to None.
        offset (int): The offset for the label position. Defaults to 15.

    Methods:
        to_plotly_annotation() -> Dict[str, Any]:
            Convert the label to Plotly annotation format.

        get_label_type() -> LabelType:
            Return the type of label. This method must be implemented by subclasses.
    """

    def __init__(
        self, position: Point3D, text: str, color: Optional[str] = None, offset: int = 15
    ) -> None:
        """Initialize the BaseLabel.

        Args:
            position (Point3D): The 3D position of the label.
            text (str): The text content of the label.
            color (Optional[str], optional): The color of the label text. Defaults to None.
            offset (int, optional): The offset for the label position. Defaults to 15.
        """
        self.position = position
        self.text = text
        self.color = color
        self.offset = offset

    def to_plotly_annotation(self) -> Dict[str, Any]:
        """Convert label to Plotly annotation format."""
        return {
            "text": self.text,
            "x": self.position.x,
            "y": self.position.y,
            "z": self.position.z,
            "showarrow": False,
            "yshift": self.offset,
            **({"font": {"color": self.color}} if self.color else {}),
        }

    @abstractmethod
    def get_label_type(self) -> LabelType:
        """Return the type of label."""
        pass


class AtomLabel(BaseLabel):
    """Class representing atom labels in molecular structure."""

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
        """Initialize the AtomLabel.

        Args:
            position (Point3D): The 3D position of the label.
            text (str): The text content of the label.
            atom_index (int): The index of the atom.
            element (str): The chemical element symbol.
            label_type (LabelType): The type of atom label.
            color (Optional[str], optional): The color of the label text. Defaults to None.
            offset (int, optional): The offset for the label position. Defaults to 15.
        """
        super().__init__(position, text, color, offset)
        self.atom_index = atom_index
        self.element = element
        self._label_type = label_type

    def get_label_type(self) -> LabelType:
        """Return the type of label for this atom."""
        return self._label_type

    @classmethod
    def from_atom(
        cls: Type["AtomLabel"],
        atom_index: int,
        element: str,
        position: Point3D,
        label_type: LabelType,
        offset: int = 15,
    ) -> "AtomLabel":
        """Create an atom label from atom data."""
        text = str(atom_index) if label_type == LabelType.ATOM_INDEX else element
        return cls(
            position=position,
            text=text,
            atom_index=atom_index,
            element=element,
            label_type=label_type,
            offset=offset,
        )


class BondLabel(BaseLabel):
    """Class representing bond labels in molecular structure."""

    def __init__(
        self,
        position: Point3D,
        bond_length: float,
        atom1_index: int,
        atom2_index: int,
        color: str = "steelblue",
        offset: int = 15,
    ) -> None:
        """Initialize the BondLabel.

        Args:
            position (Point3D): The 3D position of the label.
            bond_length (float): The length of the bond.
            atom1_index (int): The index of the first atom.
            atom2_index (int): The index of the second atom.
            color (str, optional): The color of the label text. Defaults to "steelblue".
            offset (int, optional): The offset for the label position. Defaults to 15.
        """
        super().__init__(position, f"{bond_length:.2f}", color, offset)
        self.atom1_index = atom1_index
        self.atom2_index = atom2_index
        self.bond_length = bond_length

    def get_label_type(self) -> LabelType:
        """Return the type of label for this bond."""
        return LabelType.BOND_LENGTH

    @classmethod
    def from_bond(
        cls: Type["BondLabel"],
        atom1_index: int,
        atom2_index: int,
        atom1_pos: Point3D,
        atom2_pos: Point3D,
        bond_length: float,
        color: str = "steelblue",
        offset: int = 15,
    ) -> "BondLabel":
        """Create a bond label from bond data."""
        position = Point3D.midpoint(atom1_pos, atom2_pos)
        return cls(
            position=position,
            bond_length=bond_length,
            atom1_index=atom1_index,
            atom2_index=atom2_index,
            color=color,
            offset=offset,
        )


class PlotlyButtonFactory:
    """Factory class for creating Plotly button configurations."""

    @staticmethod
    def create_toggle_button(
        label: str, visible_indices: List[int], all_indices: List[int]
    ) -> Dict[str, Any]:
        """Create a single toggle button configuration."""
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
    """Manager class for molecular structure labels with improved organization."""

    def __init__(self) -> None:
        """Initialize an empty MolecularLabelManager."""
        self._labels: List[BaseLabel] = []

    def add_label(self, label: BaseLabel) -> None:
        """Add a new label."""
        self._labels.append(label)

    def clear(self) -> None:
        """Remove all labels."""
        self._labels.clear()

    def get_labels_by_type(self, label_type: LabelType) -> List[BaseLabel]:
        """Get labels of a specific type."""
        return [label for label in self._labels if label.get_label_type() == label_type]

    def to_plotly_annotations(self) -> List[Dict[str, Any]]:
        """Convert all labels to Plotly annotation format."""
        return [label.to_plotly_annotation() for label in self._labels]

    def _get_label_indices_by_type(self, label_type: LabelType) -> List[int]:
        """Get indices of labels of a specific type."""
        return [i for i, label in enumerate(self._labels) if label.get_label_type() == label_type]

    def create_toggle_buttons(self) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
        """Create Plotly buttons for toggling label visibility in groups."""
        element_indices = self._get_label_indices_by_type(LabelType.ATOM_ELEMENT)
        index_indices = self._get_label_indices_by_type(LabelType.ATOM_INDEX)
        bond_indices = self._get_label_indices_by_type(LabelType.BOND_LENGTH)

        buttons = []
        annotations = []

        # Configuration constants
        BUTTON_LAYERS = {"distances": 1.12, "atoms": 1.065}

        # Create Bond Lengths buttons
        if bond_indices:
            buttons.append(
                {
                    "buttons": [
                        PlotlyButtonFactory.create_toggle_button(
                            "Show" if idx == 0 else "Hide",
                            bond_indices if idx == 0 else [],
                            bond_indices,
                        )
                        for idx in range(2)
                    ],
                    "direction": "right",
                    "pad": {"r": 10, "t": 10},
                    "showactive": True,
                    "active": 1,
                    "x": 0.07,
                    "y": BUTTON_LAYERS["distances"],
                    "xanchor": "left",
                    "yanchor": "top",
                    "type": "buttons",
                }
            )

            annotations.append(
                {
                    "text": "Bond Lengths",
                    "x": 0,
                    "xref": "paper",
                    "y": 1.1,
                    "yref": "paper",
                    "align": "left",
                    "showarrow": False,
                }
            )

        # Create Atom Labels buttons
        if element_indices and index_indices:
            all_atom_indices = element_indices + index_indices
            buttons.append(
                {
                    "buttons": [
                        PlotlyButtonFactory.create_toggle_button(
                            "Element", element_indices, all_atom_indices
                        ),
                        PlotlyButtonFactory.create_toggle_button(
                            "ID", index_indices, all_atom_indices
                        ),
                        PlotlyButtonFactory.create_toggle_button("Hide", [], all_atom_indices),
                    ],
                    "direction": "right",
                    "pad": {"r": 10, "t": 10},
                    "showactive": True,
                    "x": 0.07,
                    "y": BUTTON_LAYERS["atoms"],
                    "xanchor": "left",
                    "yanchor": "top",
                    "type": "buttons",
                }
            )

            annotations.append(
                {
                    "text": "Atom Labels",
                    "x": 0,
                    "xref": "paper",
                    "y": 1.045,
                    "yref": "paper",
                    "align": "left",
                    "showarrow": False,
                }
            )

        return buttons, annotations

    def update_figure_menu(self, fig: go.Figure) -> None:
        """Update figure layout with grouped toggle buttons and initial annotations."""
        annotations = self.to_plotly_annotations()

        # Set initial visibility
        for i, label in enumerate(self._labels):
            annotations[i]["visible"] = (
                isinstance(label, AtomLabel) and label.get_label_type() == LabelType.ATOM_ELEMENT
            )

        fig.update_layout(scene=dict(annotations=annotations))

        # Add button groups and their annotations
        button_groups, group_annotations = self.create_toggle_buttons()
        if button_groups:
            fig.update_layout(updatemenus=button_groups, annotations=group_annotations)
