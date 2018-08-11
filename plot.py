import plotly.graph_objs as go
from plotly.offline import offline

cpk_colors = dict(Ar='cyan', B='salmon', Ba='darkgreen', Be='darkgreen', Br='darkred', C='black', Ca='darkgreen',
                  Cl='green', Cs='violet', F='green', Fe='darkorange', Fr='violet', H='white', He='cyan',
                  I='darkviolet', K='violet', Kr='cyan', Li='violet', Mg='darkgreen', N='blue', Na='violet', Ne='cyan',
                  O='red', P='orange', Ra='darkgreen', Rb='violet', S='yellow', Sr='darkgreen', Ti='gray', Xe='cyan')
cpk_color_rest = 'pink'


def plot(adjacency_list: dict, elements: list, x_coordinates: list,
         y_coordinates: list, z_coordinates: list, plot_name: str = 'plot') -> None:
    """Creates a 3D scatter plot"""
    def atom_trace():
        """Creates an atom trace for the plot"""
        colors = [cpk_colors.get(element, cpk_color_rest) for element in elements]
        markers = dict(color=colors, line=dict(color='lightgray', width=2), size=7, symbol='circle', opacity=0.8)
        trace = go.Scatter3d(x=x_coordinates, y=y_coordinates, z=z_coordinates, mode='markers', marker=markers,
                             text=elements)
        return trace

    def bond_trace():
        """"Creates a bond trace for the plot"""
        trace = go.Scatter3d(x=[], y=[], z=[], hoverinfo='none', mode='lines',
                             marker=dict(color='grey', size=7, opacity=1))
        adjascent_atoms = ((atom, neighbour) for atom, neighbours in adjacency_list.items()
                           for neighbour in neighbours)
        for i, j in adjascent_atoms:
            trace['x'] += (x_coordinates[i], x_coordinates[j], None)
            trace['y'] += (y_coordinates[i], y_coordinates[j], None)
            trace['z'] += (z_coordinates[i], z_coordinates[j], None)
        return trace

    atoms = zip(elements, x_coordinates, y_coordinates, z_coordinates)
    annotations = [dict(text=element, x=x, y=y, z=z, showarrow=False, yshift=15) for element, x, y, z in atoms]
    data = [atom_trace(), bond_trace()]
    axis_params = dict(showgrid=False, showticklabels=False, zeroline=False, titlefont=dict(color='white'))
    layout = go.Layout(scene=dict(xaxis=axis_params, yaxis=axis_params, zaxis=axis_params, annotations=annotations),
                       margin=dict(r=0, l=0, b=0, t=0), showlegend=False)
    fig = go.Figure(data=data, layout=layout)
    offline.plot(fig, show_link=False, filename=plot_name + '.html')