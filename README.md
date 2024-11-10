# xyz2graph

[![PyPI version](https://img.shields.io/pypi/v/xyz2graph.svg)](https://pypi.org/project/xyz2graph/)
[![Python Version](https://img.shields.io/pypi/pyversions/xyz2graph.svg)](https://pypi.org/project/xyz2graph/)
[![License](https://img.shields.io/github/license/zotko/xyz2graph.svg)](https://github.com/zotko/xyz2graph/blob/master/LICENSE)

[![PyPI Downloads](https://static.pepy.tech/badge/xyz2graph)](https://pepy.tech/projects/xyz2graph)
[![GitHub Stars](https://img.shields.io/github/stars/zotko/xyz2graph)](https://github.com/zotko/xyz2graph/stargazers)
[![GitHub Forks](https://img.shields.io/github/forks/zotko/xyz2graph)](https://github.com/zotko/xyz2graph/network/members)

[!["Buy Me A Coffee"](https://img.shields.io/badge/Buy%20Me%20a%20Coffee-ffdd00?style=flat&logo=buy-me-a-coffee&logoColor=black)](https://www.buymeacoffee.com/mykola_zotko)
[![Stand With Ukraine](https://img.shields.io/badge/Stand%20With-Ukraine-FFD500?style=flat&labelColor=005BBB)](https://stand-with-ukraine.pp.ua)

A lightweight Python package for reading XYZ files and creating interactive molecular visualizations. Convert molecular structures into 3D visualizations using Plotly or analyze them as NetworkX graphs.

<div align="center">
 <img src="https://raw.githubusercontent.com/zotko/xyz2graph/main/.github/images/mol.gif" width="1024">
</div>

## Features

- Read and parse XYZ molecular structure files
- Generate interactive 3D molecular visualizations using Plotly
- Convert molecular structures to NetworkX graphs for analysis
- Command-line interface for quick visualizations

## Installation

```bash
pip install xyz2graph
```


## Requirements

- Python 3.8+
- NumPy
- Plotly
- NetworkX

## Usage

### Basic Usage

```python
from xyz2graph import MolGraph

# Create molecular graph and read XYZ file
mg = MolGraph()
mg.read_xyz('molecule.xyz')

# Generate interactive 3D visualization
fig = mg.to_plotly()
fig.show()

# Convert to NetworkX graph
G = mg.to_networkx()
```

### Command Line

```bash
# Save visualization as HTML
xyz2graph molecule.xyz

# Specify output file
xyz2graph molecule.xyz --output visualization.html

# Open directly in browser
xyz2graph molecule.xyz --browser
```

## Examples

Example XYZ files can be found in the [`examples`](examples/) directory.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation
If you use xyz2graph in your research, please cite:

```bibtex
@misc{xyz2graph,
  author       = {Zotko, Mykola},
  title        = {xyz2graph: Molecular Structure Visualization},
  year         = {2018},
  publisher    = {GitHub},
  url          = {https://github.com/zotko/xyz2graph}
}
```
