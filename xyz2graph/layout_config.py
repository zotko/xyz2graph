"""Layout configuration for molecular visualization components.

This module provides shared layout settings for molecular visualization components
including atom traces, bond traces, atom labels, bond length labels, and background controls.
"""

from typing import TypedDict


class _CommonLayout(TypedDict):
    """Common layout configuration for all buttons and annotations."""

    button_x: float
    text_x: float


class _ComponentLayout(TypedDict):
    """Layout configuration for individual component groups."""

    button_y: float
    text_y: float
    title: str


class _LayoutConfig(TypedDict):
    """Complete layout configuration for all molecular visualization components."""

    atom_traces: _ComponentLayout
    bond_traces: _ComponentLayout
    atom_labels: _ComponentLayout
    bond_lengths: _ComponentLayout
    background: _ComponentLayout
    common: _CommonLayout


LAYOUT: _LayoutConfig = {
    # Traces come first (top)
    "atom_traces": {"button_y": 1.175, "text_y": 1.155, "title": "Atoms"},
    "bond_traces": {"button_y": 1.12, "text_y": 1.1, "title": "Bonds"},
    # Then labels
    "atom_labels": {"button_y": 1.065, "text_y": 1.045, "title": "Atom Labels"},
    "bond_lengths": {"button_y": 1.01, "text_y": 0.99, "title": "Bond Lengths"},
    # Background controls
    "background": {"button_y": 0.955, "text_y": 0.935, "title": "Background"},
    # Common settings
    "common": {"button_x": 0.07, "text_x": 0},
}
