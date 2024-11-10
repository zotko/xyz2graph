# xyz2graph

[![PyPI version](https://img.shields.io/pypi/v/xyz2graph.svg)](https://pypi.org/project/xyz2graph/)
[![Python Version](https://img.shields.io/pypi/pyversions/xyz2graph.svg)](https://pypi.org/project/xyz2graph/)
[![License](https://img.shields.io/github/license/zotko/xyz2graph.svg)](https://github.com/zotko/xyz2graph/blob/master/LICENSE)

[![PyPI Downloads](https://static.pepy.tech/badge/xyz2graph)](https://pepy.tech/projects/xyz2graph)
[![GitHub Stars](https://img.shields.io/github/stars/zotko/xyz2graph)](https://github.com/zotko/xyz2graph/stargazers)
[![GitHub Forks](https://img.shields.io/github/forks/zotko/xyz2graph)](https://github.com/zotko/xyz2graph/network/members)

[!["Buy Me A Coffee"](https://img.shields.io/badge/Buy%20Me%20a%20Coffee-ffdd00?style=flat&logo=buy-me-a-coffee&logoColor=black)](https://www.buymeacoffee.com/mykola_zotko)
[![Stand With Ukraine](https://img.shields.io/badge/Stand%20With-Ukraine-FFD500?style=flat&labelColor=005BBB)](https://stand-with-ukraine.pp.ua)

Welcome to xyz2graph's documentation! This Python package provides tools for reading XYZ molecular structure files and creating interactive 3D visualizations.

<div align="center">
 <img src="https://raw.githubusercontent.com/zotko/xyz2graph/main/.github/images/mol.gif" width="1024">
</div>

## Key Features
- Read and parse XYZ molecular structure files
- Generate interactive 3D molecular visualizations using Plotly
- Convert molecular structures to NetworkX graphs
- Command-line interface for quick visualizations

## Requirements

- Python 3.8+
- NumPy
- Plotly
- NetworkX

## Installation
```bash
pip install xyz2graph
```

## Quick Example
```python
from xyz2graph import MolGraph

mg = MolGraph()
mg.read_xyz('molecule.xyz')
fig = mg.to_plotly()
fig.show()
```

## Quick Navigation
- [Getting Started](getting-started.md)
- API Documentation
    - [MolGraph](api/molgraph.md)
    - [Command Line](api/cli.md)

## Contributing
Contributions are welcome! Please feel free to submit a Pull Request.

## License
This project is licensed under the MIT License - see the [LICENSE](https://github.com/zotko/xyz2graph/blob/master/LICENSE) file for details.
