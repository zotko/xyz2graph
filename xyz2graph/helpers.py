import plotly.graph_objs as go
import networkx as nx

from xyz2graph import MolGraph

cpk_colors = dict(Ar='cyan', B='salmon', Ba='darkgreen', Be='darkgreen', Br='darkred', C='black', Ca='darkgreen',
                  Cl='green', Cs='violet', F='green', Fe='darkorange', Fr='violet', H='white', He='cyan',
                  I='darkviolet', K='violet', Kr='cyan', Li='violet', Mg='darkgreen', N='blue', Na='violet', Ne='cyan',
                  O='red', P='orange', Ra='darkgreen', Rb='violet', S='yellow', Sr='darkgreen', Ti='gray', Xe='cyan')
cpk_color_rest = 'pink'


def to_plotly_figure(graph: MolGraph) -> go.Figure:
    """Creates a Plotly figure."""

    def atom_trace():
        """Creates an atom trace for the plot."""
        colors = [cpk_colors.get(element, cpk_color_rest) for element in graph.elements]
        markers = dict(color=colors, line=dict(color='lightgray', width=2), size=7, symbol='circle', opacity=0.8)
        trace = go.Scatter3d(x=graph.x, y=graph.y, z=graph.z, mode='markers',
                             marker=markers,
                             text=graph.elements)
        return trace

    def bond_trace():
        """"Creates a bond trace for the plot."""
        trace = go.Scatter3d(x=[], y=[], z=[], hoverinfo='none', mode='lines',
                             marker=dict(color='grey', size=7, opacity=1))
        adjascent_atoms = ((atom, neighbour) for atom, neighbours in graph.adj_list.items()
                           for neighbour in neighbours)
        for i, j in adjascent_atoms:
            trace['x'] += (graph.x[i], graph.x[j], None)
            trace['y'] += (graph.y[i], graph.y[j], None)
            trace['z'] += (graph.z[i], graph.z[j], None)
        return trace

    annotations_elements = [dict(text=element, x=x, y=y, z=z, showarrow=False, yshift=15)
                            for element, (x, y, z) in graph]
    annotations_indices = [dict(text=number, x=x, y=y, z=z, showarrow=False, yshift=15)
                           for number, (_, (x, y, z)) in enumerate(graph)]

    annotations_bonds = []
    for (i, j), length in graph.bond_lengths.items():
        x = (graph.x[i] + graph.x[j]) / 2
        y = (graph.y[i] + graph.y[j]) / 2
        z = (graph.z[i] + graph.z[j]) / 2
        annotations_bonds.append(dict(text=round(length, 2), x=x, y=y, z=z,
                                      showarrow=False, yshift=15, font=dict(color="steelblue")))

    updatemenus = list([
        dict(buttons=list([
            dict(label=' Elements',
                 method='relayout',
                 args=[{'scene.annotations': annotations_elements}]),
            dict(label=' Elements & Bond Lengths',
                 method='relayout',
                 args=[{'scene.annotations': annotations_elements + annotations_bonds}]),
            dict(label='Indices',
                 method='relayout',
                 args=[{'scene.annotations': annotations_indices}]),
            dict(label='Indices & Bond Lengths',
                 method='relayout',
                 args=[{'scene.annotations': annotations_indices + annotations_bonds}]),
            dict(label='Bond Lengths',
                 method='relayout',
                 args=[{'scene.annotations': annotations_bonds}]),
            dict(label='Hide All',
                 method='relayout',
                 args=[{'scene.annotations': []}])
        ]),
            direction='down',
            xanchor='left',
            yanchor='top'
        ),
    ])

    data = [atom_trace(), bond_trace()]
    axis_params = dict(showgrid=False, showbackground=False, showticklabels=False, zeroline=False,
                       titlefont=dict(color='white'))
    layout = dict(scene=dict(xaxis=axis_params, yaxis=axis_params, zaxis=axis_params,
                             annotations=annotations_elements),
                  margin=dict(r=0, l=0, b=0, t=0), showlegend=False, updatemenus=updatemenus)
    figure = go.Figure(data=data, layout=layout)

    return figure


def to_networkx_graph(graph: MolGraph) -> nx.Graph:
    """Creates a NetworkX graph.
    Atomic elements and coordinates are added to the graph as node attributes 'element' and 'xyz" respectively.
    Bond lengths are added to the graph as edge attribute 'length''"""
    G = nx.Graph(graph.adj_list)
    node_attrs = {num: {'element': element, 'xyz': xyz} for num, (element, xyz) in enumerate(graph)}
    nx.set_node_attributes(G, node_attrs)
    edge_attrs = {edge: {'length': length} for edge, length in graph.bond_lengths.items()}
    nx.set_edge_attributes(G, edge_attrs)
    return G
