# Getting Started

This guide will help you get started with xyz2graph.

## Installation

Install xyz2graph using pip:

```bash
pip install xyz2graph
```

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

### Using the Command Line

For quick visualizations, you can use the command-line interface:

```bash
xyz2graph molecule.xyz
```

This will create an HTML file with the visualization in your current directory.

## Customization

### Changing Colors

You can customize the appearance of your molecular visualization by changing element colors:

```python
mg = MolGraph()
mg.read_xyz('molecule.xyz')

# Change specific element colors
mg.set_element_color('O', 'pink')

# Set default color for unspecified elements
mg.set_default_color('gray')

# Create visualization with custom colors
fig = mg.to_plotly()
fig.show()
```

### Adjusting Atomic Radii

Modify atomic radii to adjust the visual representation:

```python
# Change radius for specific elements
mg.set_element_radius('C', 0.75)
mg.set_element_radius('H', 0.25)
```

## Network Analysis

Convert your molecular structure to a NetworkX graph for analysis:

```python
# Convert to NetworkX graph
G = mg.to_networkx()

# Access molecular properties
print(f"Number of atoms: {len(G.nodes)}")
print(f"Number of bonds: {len(G.edges)}")
```

## XYZ File Format

The XYZ file format is a simple text format for molecular structures:

```
3
Water
O  0.000000  0.000000  0.000000
H  0.758602  0.000000  0.504284
H  0.758602  0.000000 -0.504284
```

Where:

- First line: Number of atoms
- Second line: Comment or title
- Following lines: Element and coordinates (x, y, z) in Angstroms

## Next Steps

- See the [API documentation](api/molgraph.md) for detailed reference
- Learn about [command-line options](api/cli.md) for batch processing

## Common Issues

### File Not Found

If you get a `FileNotFoundError`, check that:

- The XYZ file exists in the specified path
- You have read permissions for the file
- The path is correct relative to your working directory

### Invalid File Format

If you get a `ValueError` when reading an XYZ file:

- Verify the file follows the XYZ format
- Check for missing or extra whitespace
- Ensure coordinates are valid numbers
- Verify element symbols are valid
- If your file doesn't have the standard header (atom count and comment),
use the `xyz_start` parameter to specify where coordinates begin:
  ```python
  # Start reading from first line (no header)
  mg.read_xyz("molecule.xyz", xyz_start=0)

  # Start reading from line 3
  mg.read_xyz("molecule.xyz", xyz_start=3)
  ```

For more specific issues raise an issue on our [GitHub repository](https://github.com/zotko/xyz2graph/issues).
