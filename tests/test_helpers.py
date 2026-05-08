"""Tests for the deprecated helpers module."""

import warnings
from typing import Any, Callable

import networkx as nx
import plotly.graph_objects as go
import pytest
from xyz2graph.graph import MolGraph
from xyz2graph.helpers import to_networkx_graph, to_plotly_figure


@pytest.mark.parametrize(
    ("func", "expected_type"),
    [
        (to_networkx_graph, nx.Graph),
        (to_plotly_figure, go.Figure),
    ],
)
def test_deprecated_helpers_emit_warning(
    water_molecule: MolGraph, func: Callable[[MolGraph], Any], expected_type: type
) -> None:
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = func(water_molecule)

    assert isinstance(result, expected_type)
    assert any(issubclass(w.category, DeprecationWarning) for w in caught)
