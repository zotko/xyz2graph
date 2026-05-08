"""Tests for the visualization module."""

from typing import cast

import plotly.graph_objects as go
import pytest
from xyz2graph.graph import MolGraph
from xyz2graph.molecule import Atom
from xyz2graph.visualization import (
    DEFAULT_CONFIG,
    AtomShape,
    BondShape,
    ConfigurationManager,
    SceneConfig,
    StyleConfig,
    VisualizationConfig,
    VisualizationManager,
    WatermarkConfig,
    create_visualization,
)


@pytest.fixture
def lone_atom() -> MolGraph:
    mol = MolGraph()
    mol.atoms = [Atom("He", 0.0, 0.0, 0.0, 0, radius=0.28)]
    return mol


@pytest.mark.parametrize(
    ("label_type", "expected"),
    [
        ("indices", {"0", "1", "2"}),
        ("none", {""}),
    ],
)
def test_atom_shape_label_type(
    water_molecule: MolGraph, label_type: str, expected: set[str]
) -> None:
    shape = AtomShape(atoms=water_molecule.atoms, colors=water_molecule.cpk_colors)
    traces = shape.to_trace(
        config=DEFAULT_CONFIG["style"]["atoms"],
        coordinate_round_digits=3,
        label_type=label_type,
    )
    rendered = {t for trace in traces for t in trace.text}
    assert rendered == expected


def test_atom_shape_unknown_element_uses_pink_default() -> None:
    atoms = [Atom("X", 0.0, 0.0, 0.0, 0)]
    shape = AtomShape(atoms=atoms, colors={})
    traces = shape.to_trace(
        config=DEFAULT_CONFIG["style"]["atoms"],
        coordinate_round_digits=3,
    )
    assert traces[0].marker.color == "pink"


def test_bond_shape_to_trace_emits_line_and_text_pairs(water_molecule: MolGraph) -> None:
    shape = BondShape(bonds=water_molecule.bonds)
    traces = shape.to_trace(
        config=DEFAULT_CONFIG["style"]["bonds"],
        coordinate_round_digits=3,
        markdown_round_digits=2,
    )
    assert len(traces) == 2
    assert traces[0].mode == "lines"
    assert traces[1].mode == "text"


def test_configuration_manager_defaults() -> None:
    cm = ConfigurationManager()
    assert cm.style == DEFAULT_CONFIG["style"]
    assert cm.scene == DEFAULT_CONFIG["scene"]
    assert cm.buttons == DEFAULT_CONFIG["buttons"]
    assert cm.atoms == DEFAULT_CONFIG["style"]["atoms"]
    assert cm.bonds == DEFAULT_CONFIG["style"]["bonds"]
    assert cm.watermark == DEFAULT_CONFIG["style"]["watermark"]


def test_to_plotly_full_figure(water_molecule: MolGraph) -> None:
    fig = water_molecule.to_plotly()
    assert isinstance(fig, go.Figure)
    assert len(fig.data) == 6
    assert len(fig.layout.updatemenus) == 3


def test_create_visualization_no_bonds_skips_bond_controls(lone_atom: MolGraph) -> None:
    fig = create_visualization(lone_atom)
    assert isinstance(fig, go.Figure)
    assert len(fig.layout.updatemenus) == 2


def test_create_visualization_with_long_comment_is_shortened(water_molecule: MolGraph) -> None:
    water_molecule.comment = "x" * 500
    fig = create_visualization(water_molecule)
    annotation_texts = [a.text for a in fig.layout.annotations]
    truncated = next(t for t in annotation_texts if t and t.endswith("..."))
    assert len(truncated) <= DEFAULT_CONFIG["style"]["information"]["comment"]["max_length"]


def test_create_visualization_axis_visible_off(water_molecule: MolGraph) -> None:
    scene = cast(SceneConfig, {**DEFAULT_CONFIG["scene"], "showbackground": False})
    config: VisualizationConfig = {
        "style": DEFAULT_CONFIG["style"],
        "scene": scene,
        "buttons": DEFAULT_CONFIG["buttons"],
    }
    manager = VisualizationManager(water_molecule, config)
    axis_config = manager._create_axis_config()
    assert axis_config["xaxis"]["gridcolor"] == "white"
    assert axis_config["xaxis"]["backgroundcolor"] is None


def test_create_visualization_no_watermark(water_molecule: MolGraph) -> None:
    watermark = cast(WatermarkConfig, {**DEFAULT_CONFIG["style"]["watermark"], "show": False})
    style = cast(StyleConfig, {**DEFAULT_CONFIG["style"], "watermark": watermark})
    config: VisualizationConfig = {
        "style": style,
        "scene": DEFAULT_CONFIG["scene"],
        "buttons": DEFAULT_CONFIG["buttons"],
    }
    fig = create_visualization(water_molecule, config)
    assert all("xyz2graph" not in (a.text or "") for a in fig.layout.annotations)


def test_create_visualization_empty_comment_skipped(water_molecule: MolGraph) -> None:
    water_molecule.comment = ""
    fig = create_visualization(water_molecule)
    info_texts = [
        a.text for a in fig.layout.annotations if a.text in {"H2O"} or "xyz2graph" in (a.text or "")
    ]
    assert "H2O" in info_texts
