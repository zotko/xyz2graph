# xyz2graph

[**xyz2graph**](https://github.com/zotko/xyz2graph) is a Python package for reading of .xyz files and constructing of molecular graphs from atomic coordinates. The molecular graph can be converted into [NetworkX](https://networkx.github.io) graph or [Plotly](https://plot.ly) figure for 3D visualization in a browser window or in a [Jupyter notebook](https://jupyter.org).

<p align="center">
  <img src="picture.png",  width="512">
</p>

## Requirements
 * **Python 3**
 * [NetworkX](https://networkx.github.io)
 * [Plotly](https://plot.ly)
 
## Installation
`!pip install -i https://test.pypi.org/simple/ xyz2graph`

## Usage
```
from xyz2graph import MolGraph, to_networkx_graph, to_plotly_figure
from plotly.offline import offline

# Create the MolGraph object
mg = MolGraph()

# Read the data from the .xyz file
mg.read_xyz('file_path')

# Create the Plotly figure object
fig = to_plotly_figure(mg)

# Plot the figure
offline.plot(fig)

# Convert the molecular graph to the NetworkX graph
G = to_networkx_graph(mg)
```

## Usage in Jupyter Notebook

```
from xyz2graph import MolGraph, to_networkx_graph, to_plotly_figure
from plotly.offline import init_notebook_mode, iplot

# Initiate the Plotly notebook mode
init_notebook_mode(connected=True)

# Create the MolGraph object
mg = MolGraph()

# Read the data from the .xyz file
mg.read_xyz('file_path')

# Create the Plotly figure object
fig = to_plotly_figure(mg)

# Plot the figure
iplot(fig)
```

### [Example of usage](https://www.kaggle.com/mykolazotko/xyz2graph-xyz-file-to-molecular-graph?scriptVersionId=20042579)
