import warnings

import networkx as nx
import plotly.graph_objs as go

from .graph import MolGraph


def to_networkx_graph(graph: MolGraph) -> nx.Graph:
    """Deprecated: Use MolGraph.to_networkx() method instead."""
    warnings.warn(
        "to_networkx_graph is deprecated and will be removed in version 4.0.0. "
        + "Use MolGraph.to_networkx() instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return graph.to_networkx()


def to_plotly_figure(graph: MolGraph) -> go.Figure:
    """Deprecated: Use MolGraph.to_plotly() method instead."""
    warnings.warn(
        "to_plotly_figure is deprecated and will be removed in version 4.0.0. "
        + "Use MolGraph.to_plotly() instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return graph.to_plotly()
