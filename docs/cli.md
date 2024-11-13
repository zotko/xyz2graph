# Command Line Interface

## Usage

```bash
xyz2graph [-h] [-o OUTPUT] [-b] [--debug] [-f FILTER] xyz_file
```

## Options

**xyz_file**  
Input XYZ file path containing molecular coordinates (required).

**-o, --output**  
Custom output path for the HTML visualization. Defaults to input filename with .html extension.

**-b, --browser**  
Launch visualization directly in default web browser instead of saving to file.

**-f, --filter**  
Filter out specific atoms using comma-separated indices or element symbols (e.g., "1,2,H,O,5" removes atoms
at positions 1,2,5 and all H,O atoms). Element symbols must be properly capitalized.

**--debug**  
Enable verbose logging for troubleshooting and development.

**-h, --help**  
Display command-line options and usage information.

## Examples

```bash
# Basic usage - saves visualization as HTML in the same directory
xyz2graph molecule.xyz

# Save visualization with a specific name
xyz2graph molecule.xyz --output my-molecule.html

# Open directly in browser without saving
xyz2graph molecule.xyz --browser

# Remove hydrogens and oxygens from visualization
xyz2graph molecule.xyz --filter "H,O"

# Remove specific atoms by index and element
xyz2graph molecule.xyz --filter "1,2,H,O,5"
```
