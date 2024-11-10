# MolGraph
The `MolGraph` class is the core component of xyz2graph, providing functionality to read molecular structure data and convert it to different representations.

## Basic Usage
```python
from xyz2graph import MolGraph

# Create a new molecular graph
mg = MolGraph()

# Read molecular structure from XYZ file
mg.read_xyz('molecule.xyz')

# Create visualization
fig = mg.to_plotly()
fig.show()
```

## Methods

### `read_xyz`
Reads and parses molecular structure data from an XYZ file.
```python
mg.read_xyz('molecule.xyz')
```

- `file_path` (str or Path): Path to the XYZ file
- `xyz_start` (int): Line number where coordinates start (default: 2)
- `validate` (bool): Validate atom count against file header (default: False)

### `to_plotly`
Creates a plotly figure with interactive 3D visualization of the molecule.
```python
fig = mg.to_plotly()
```

- Returns: Plotly Figure object for 3D visualization

### `to_networkx`
Creates a NetworkX graph representation of the molecule.
```python
G = mg.to_networkx()
```

- Returns: NetworkX Graph object with node and edge attributes

### `formula`
Returns the molecular formula using Hill notation convention.
```python
formula = mg.formula()
print(formula)  # e.g., "CH4O"
```

- Returns: str, molecular formula

### `set_element_radius`
Sets the reference radius for a specific element.
```python
mg.set_element_radius('C', 0.75)  # Set carbon atom radius
```

- `element` (str): Chemical element symbol
- `radius` (float): New radius value in Angstroms

### `set_element_color`
Sets the CPK color for a specific element.
```python
mg.set_element_color('O', 'red')  # Set oxygen atoms to red
```

- `element` (str): Chemical element symbol
- `color` (str): Color name or code accepted by Plotly

### `set_default_color`
Sets the default color for elements without specified colors.
```python
mg.set_default_color('gray')
```

- `color` (str): Color name or code accepted by Plotly

## Properties

### `elements`
List of chemical element symbols in the molecule

### `bond_lengths`
Dictionary mapping pairs of connected atoms to their bond lengths

### `adj_list`
Dictionary mapping atom indices to their connected neighbors

### `atomic_radii`
List of atomic radii for each atom in the molecule

### `x, y, z`
Lists of coordinates for each atom

These properties are automatically updated when reading an XYZ file or modifying the molecular structure.
