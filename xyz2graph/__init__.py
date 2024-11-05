import warnings

from .xyz2graph import MolGraph


def to_networkx_graph(graph):
    warnings.warn(
        "to_networkx_graph is deprecated and will be removed in version 1.0.0. "
        "Use MolGraph.to_networkx() instead. "
        "See: https://github.com/YOUR_USERNAME/YOUR_REPO/blob/main/MIGRATION.md",
        DeprecationWarning,
        stacklevel=2,
    )
    return graph.to_networkx()


def to_plotly_figure(graph):
    warnings.warn(
        "to_plotly_figure is deprecated and will be removed in version 1.0.0. "
        "Use MolGraph.to_plotly() instead. "
        "See: https://github.com/YOUR_USERNAME/YOUR_REPO/blob/main/MIGRATION.md",
        DeprecationWarning,
        stacklevel=2,
    )
    return graph.to_plotly()


__all__ = ["MolGraph", "to_networkx_graph", "to_plotly_figure"]
