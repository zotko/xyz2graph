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

### Advanced Reading Options

Control how files are read and validated:

```python
# Read raw coordinates (no header)
mg.read_xyz('coords.xyz', xyz_start=0)

# Read coordinates from custom start line
mg.read_xyz('custom.xyz', xyz_start=3)

# Validate coordinate count against header
mg.read_xyz('molecule.xyz', validate=True)
```

### Network Analysis

Convert your molecular structure to a NetworkX graph for analysis:

```python
# Convert to NetworkX graph
G = mg.to_networkx()

# Calculate graph properties
print(f"Graph density: {nx.density(G)}")
```

### Visualization Settings

You can customize the appearance of your molecular visualization:

```python
config = {
    "atom_size": 7,
    "atom_opacity": 0.8,
    "atom_border_color": "lightgray",
    "atom_border_width": 2,
    "bond_color": "grey",
    "bond_width": 2,
    "show_grid": False,
    "label_offset": 15,
    "bond_label_color": "steelblue"
}

# Create visualization with custom settings
fig = mg.to_plotly(config=config)
fig.show()
```

### Structure Manipulation

```python
# Remove all hydrogen atoms
mg.filter(elements=['H'], inplace=True)

# Remove multiple element types
mg.filter(elements=['H', 'O'], inplace=True)

# Remove all hydrogen atoms and atoms at indices 1 and 5
new_mg = mg.filter(elements=['H'], indices=[1,5])
```

## Using the Command Line

For quick visualizations, you can use the command-line interface:

```bash
# Basic usage - saves visualization as HTML in the same directory
xyz2graph molecule.xyz

# Save visualization to specific path
xyz2graph molecule.xyz --output viz.html

# Open directly in browser without saving
xyz2graph molecule.xyz --browser

# Remove hydrogen atoms from visualization
xyz2graph molecule.xyz --filter "H"

# Remove specific atoms by index (0-based)
xyz2graph molecule.xyz --filter "0,1,2"

# Combine element and index filtering
xyz2graph molecule.xyz --filter "H,O,1,2,5"
```

## Next Steps

- See the [API documentation](python.md) for detailed reference
- Learn about [command-line options](cli.md) for batch processing

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
  # Start reading coordinates from first line (no header)
  mg.read_xyz("molecule.xyz", xyz_start=0)

  # Start reading coordinates from line 3
  mg.read_xyz("molecule.xyz", xyz_start=3)
  ```

For more specific issues raise an issue on our [GitHub repository](https://github.com/zotko/xyz2graph/issues).
