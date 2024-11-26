# Getting Started

This guide will help you get started with xyz2graph.

## Basic Usage

### Creating a Molecular Visualization

The most common use case is to create an interactive 3D visualization of a molecular structure:

```python
from xyz2graph import MolGraph

# Create a molecular graph
mg = MolGraph()

# Read structure from XYZ file
mg.read_xyz('molecule.xyz')

# Generate and display visualization
fig = mg.to_plotly()
fig.show()
```

### Network Analysis

Convert your molecular structure to a NetworkX graph for analysis:

```python
import networkx as nx

# Convert to NetworkX graph
G = mg.to_networkx()

# Calculate graph properties
print(f"Graph density: {nx.density(G)}")
```

### Structure Manipulation

Remove atoms for better visualization.

```python
# Remove all hydrogen atoms, modifying the original structure
mg.remove(elements=['H'], inplace=True)

# Remove both hydrogen and oxygen atoms, modifying the original structure
mg.remove(elements=['H', 'O'], inplace=True)

# Create a new molecule by removing hydrogens and specific atoms
# The original molecule remains unchanged since inplace=False (default)
new_mg = mg.remove(elements=['H'], indices=[1, 5])
```

## Using the Command Line

For quick visualizations, you can use the command-line interface:

```bash
# Basic usage - saves visualization as HTML in the same directory
xyz2graph molecule.xyz

# Save visualization to specific path
xyz2graph molecule.xyz -o viz.html

# Open directly in browser without saving
xyz2graph molecule.xyz -b

# Remove hydrogen atoms from visualization
xyz2graph molecule.xyz -r "H"

# Remove specific atoms by index (0-based)
xyz2graph molecule.xyz -r "0,1,2"

# Combine element and index removal
xyz2graph molecule.xyz -r "H,O,1,2,5"
```

## Next Steps

- See the [API documentation](python.md) for detailed reference
- Learn about [command-line options](cli.md) for batch processing

For more specific issues raise an issue on our [GitHub repository](https://github.com/zotko/xyz2graph/issues).
