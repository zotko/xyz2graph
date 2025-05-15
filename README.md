# xyz2graph

[![PyPI version](https://img.shields.io/pypi/v/xyz2graph.svg)](https://pypi.org/project/xyz2graph/)
[![Python Version](https://img.shields.io/pypi/pyversions/xyz2graph.svg)](https://pypi.org/project/xyz2graph/)
[![License](https://img.shields.io/github/license/zotko/xyz2graph.svg)](https://github.com/zotko/xyz2graph/blob/master/LICENSE)
[![Documentation](https://img.shields.io/badge/docs-mkdocs-blue)](https://zotko.github.io/xyz2graph)
[![DOI](https://zenodo.org/badge/144382005.svg)](https://doi.org/10.5281/zenodo.14569337)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/zotko/xyz2graph/main?urlpath=%2Fdoc%2Ftree%2Fbinder%2Fxyz2graph.ipynb)


[![PyPI Downloads](https://static.pepy.tech/badge/xyz2graph/month)](https://pepy.tech/projects/xyz2graph)
[![GitHub Stars](https://img.shields.io/github/stars/zotko/xyz2graph)](https://github.com/zotko/xyz2graph/stargazers)
[![GitHub Forks](https://img.shields.io/github/forks/zotko/xyz2graph)](https://github.com/zotko/xyz2graph/network/members)

[![Stand With Ukraine](https://raw.githubusercontent.com/vshymanskyy/StandWithUkraine/main/badges/StandWithUkraine.svg)](https://stand-with-ukraine.pp.ua)

A Python package to convert XYZ molecular files into NetworkX graphs with interactive 3D visualization using Plotly.

Try it live 🚀: [Binder](https://mybinder.org/v2/gh/zotko/xyz2graph/main?urlpath=%2Fdoc%2Ftree%2Fbinder%2Fxyz2graph.ipynb) | [GESIS Notebooks](https://notebooks.gesis.org/binder/v2/gh/zotko/xyz2graph/main?urlpath=%2Fdoc%2Ftree%2Fbinder%2Fxyz2graph.ipynb)

## Features

- [Interactive 3D molecular visualization using Plotly](https://zotko.github.io/xyz2graph/demo)
- NetworkX graph conversion for analysis
- Command-line interface

## Installation

```bash
pip install xyz2graph
```

## Requirements

- Python 3.8+
- Dependencies: NumPy, Plotly, NetworkX

## Quick Start

```python
import numpy as np
from xyz2graph import MolGraph

# Create molecular graph and read XYZ file
mg = MolGraph()
mg.read_xyz('molecule.xyz')

# Generate interactive 3D visualization
fig = mg.to_plotly()
fig.show()

# Convert to NetworkX graph
G = mg.to_networkx()

G.nodes[0]
# Output: {'element': 'C', 'xyz': (0.1718396797, 1.4440789224, 0.2473852864)}

G.edges[(0, 1)]
# Output: {'length': np.float64(1.49623)}
```

## Command Line

```bash
# Save visualization as HTML
xyz2graph molecule.xyz

# Specify output file
xyz2graph molecule.xyz --output viz.html

# Open directly in browser
xyz2graph molecule.xyz --browser
```

## Documentation

Read the [documentation](https://zotko.github.io/xyz2graph) for guides, API reference, and examples.

## Help & Discussion

🪲 [Report a bug](https://github.com/zotko/xyz2graph/issues)  
✨ [Request a feature](https://github.com/zotko/xyz2graph/discussions)

## Contributing

Contributions are welcome! Please see the [Contributing Guide](https://github.com/zotko/xyz2graph/tree/main/CONTRIBUTING.md) for guidelines.
