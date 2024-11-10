# Command Line Interface

The package provides a command-line interface for quick visualization of XYZ files.

## Basic Usage
```bash
# Basic usage - saves visualization as HTML in the same directory
xyz2graph molecule.xyz

# Specify custom output file
xyz2graph molecule.xyz --output visualization.html

# Open directly in browser without saving
xyz2graph molecule.xyz --browser

# Enable debug logging
xyz2graph molecule.xyz --debug
```

## Options
- `xyz_file`: Path to the input XYZ file (required)
- `-o, --output`: Output HTML file path (default: based on input filename)
- `-b, --browser`: Open visualization in browser without saving
- `--debug`: Enable debug logging
