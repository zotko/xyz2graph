# xyz2graph

[![PyPI version](https://img.shields.io/pypi/v/xyz2graph.svg)](https://pypi.org/project/xyz2graph/)
[![Python Version](https://img.shields.io/pypi/pyversions/xyz2graph.svg)](https://pypi.org/project/xyz2graph/)
[![License](https://img.shields.io/github/license/zotko/xyz2graph.svg)](https://github.com/zotko/xyz2graph/blob/master/LICENSE)
[![Documentation](https://img.shields.io/badge/docs-mkdocs-blue)](https://zotko.github.io/xyz2graph)
[![DOI](https://zenodo.org/badge/144382005.svg)](https://doi.org/10.5281/zenodo.14569337)

[![PyPI Downloads](https://static.pepy.tech/badge/xyz2graph/month)](https://pepy.tech/projects/xyz2graph)
[![GitHub Stars](https://img.shields.io/github/stars/zotko/xyz2graph)](https://github.com/zotko/xyz2graph/stargazers)
[![GitHub Forks](https://img.shields.io/github/forks/zotko/xyz2graph)](https://github.com/zotko/xyz2graph/network/members)

[![Stand With Ukraine](https://raw.githubusercontent.com/vshymanskyy/StandWithUkraine/main/badges/StandWithUkraine.svg)](https://stand-with-ukraine.pp.ua)

A Python package to convert XYZ molecular files into NetworkX graphs with interactive 3D visualization using Plotly.

<a href="https://zotko.github.io/xyz2graph/demo" target="_blank">Try it live 🚀</a>

## Features

- Interactive 3D molecular visualization using Plotly
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
from xyz2graph import MolGraph

# Create molecular graph and read XYZ file
mg = MolGraph()
mg.read_xyz('molecule.xyz')

# Convert to NetworkX graph
G = mg.to_networkx()

# Generate interactive 3D visualization
fig = mg.to_plotly()
fig.show()
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

## Citation

If you use xyz2graph in your research, please cite:

```bibtex
@misc{zotko2018xyz2graph,
  author       = {Zotko, Mykola},
  title        = {xyz2graph: Molecular Structure Visualization},
  year         = {2018},
  publisher    = {GitHub},
  url          = {https://github.com/zotko/xyz2graph}
}
```

<p align="center">
  <a href="https://www.buymeacoffee.com/mykola_zotko">
    <img src="https://www.buymeacoffee.com/assets/img/custom_images/yellow_img.png" alt="Buy Me A Coffee">
  </a>
</p>
