#!/usr/bin/env python3
"""Command-line interface for generating 3D molecular visualizations from XYZ files."""

import argparse
import logging
import sys
import tempfile
from pathlib import Path
from typing import List, Tuple, Union

import plotly.offline as offline
from plotly.io import write_html

from .graph import MolGraph
from .logging import logger


def parse_remove_arg(remove_str: str) -> Tuple[List[int], List[str]]:
    """Parse remove argument string into indices and elements.

    Args:
        remove_str: Comma-separated values for indices and elements
            (e.g., "1,2,H,O" removes atoms at indices 1,2 and all H,O atoms)

    Returns:
        Tuple of (indices_list, elements_list)

    Example:
        "1,2,H,O,5" -> ([1, 2, 5], ["H", "O"])
    """
    indices = []
    elements = []

    for item in remove_str.split(","):
        item = item.strip()
        if not item:
            continue
        try:
            # Try to convert to integer
            indices.append(int(item))
        except ValueError:
            # If not an integer, treat as element symbol
            if item.isalpha():
                elements.append(item)
            else:
                logger.warning(f"Ignoring invalid value: {item}")

    return indices, elements


def format_remove_message(indices: List[int], elements: List[str]) -> str:
    """Format the message showing which atoms will be removed.

    Args:
        indices: List of atom indices to remove
        elements: List of element symbols to remove

    Returns:
        str: Formatted message describing atoms to be removed
    """
    messages = []
    if indices:
        messages.append(f"Removing atoms at indices: {', '.join(map(str, indices))}")
    if elements:
        messages.append(f"Removing elements: {', '.join(elements)}")
    return "; ".join(messages)


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate 3D molecular visualizations from XYZ files"
    )
    parser.add_argument("xyz_file", type=str, help="Path to the input XYZ file")
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        help="Output HTML file path (default: based on input filename)",
    )
    parser.add_argument(
        "-b",
        "--browser",
        action="store_true",
        help="Open the visualization in browser without saving the file",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Enable debug logging",
    )
    remove_group = parser.add_mutually_exclusive_group()
    remove_group.add_argument(
        "-r",
        "--remove",
        type=str,
        help="Atoms to remove using comma-separated values for indices and/or elements. "
        "Make sure to use proper element capitalization (e.g., 'H' not 'h'). "
        "Example: '1,2,H,O,5' removes atoms at indices 1,2,5 and all H,O atoms.",
    )

    return parser.parse_args()


def generate_output_path(input_path: str, output_path: Union[str, None]) -> Path:
    """Generate the output file path based on input path."""
    if output_path:
        return Path(output_path)
    return Path(input_path).with_suffix(".html")


def log_molecule_state(mol: MolGraph, msg: str = "") -> None:
    """Log information about the molecule's current state.

    Args:
        mol: MolGraph instance to analyze
        msg: Optional message prefix

    Returns:
        None: Logs molecule formula and optionally detailed atom info if debug enabled
    """
    logger.info(f"{msg + ': ' if msg else ''}{mol.formula()}")

    if logger.level <= logging.DEBUG:
        for atom in mol.atoms:
            logger.debug(f"{atom.index}. {atom.element} ({atom.x:.3f}, {atom.y:.3f}, {atom.z:.3f})")


def main() -> None:
    """Log information about the molecule's current state.

    Args:
        mol: MolGraph instance to analyze
        msg: Optional message prefix

    Returns:
        None: Logs molecule formula and optionally detailed atom info if debug enabled
    """
    args = parse_args()

    if args.debug:
        logger.setLevel(logging.DEBUG)
        logger.debug(f"Arguments: {args}")

    # Create MolGraph instance and read XYZ file
    mg = MolGraph()
    try:
        logger.info(f"Reading XYZ file: {args.xyz_file}")
        mg.read_xyz(args.xyz_file)
        log_molecule_state(mg, "Initial molecule")

        # Handle atom removal
        if args.remove:
            indices, elements = parse_remove_arg(args.remove)
            if indices or elements:
                logger.info(format_remove_message(indices, elements))
                mg.remove(indices=indices or None, elements=elements or None, inplace=True)
                log_molecule_state(mg, "Modified molecule")

    except (FileNotFoundError, ValueError, IndexError) as e:
        logger.error(f"Error processing molecule: {e}")
        sys.exit(1)

    # Generate figure
    fig = mg.to_plotly()

    # Save/display the visualization
    try:
        if args.browser:
            logger.info("Opening visualization in browser")
            # Create a temporary file that will be automatically cleaned up
            with tempfile.NamedTemporaryFile(
                prefix=Path(args.xyz_file).stem + "_", suffix=".html", delete=False
            ) as tmp:
                offline.plot(fig, filename=tmp.name, auto_open=True)
                logger.debug(f"Created temporary file: {tmp.name}")
        else:
            output_path = generate_output_path(args.xyz_file, args.output)
            logger.info(f"Saving visualization to: {output_path}")
            write_html(
                fig,
                str(output_path),
                config={
                    "toImageButtonOptions": {
                        "filename": Path(output_path).stem,
                        "format": "svg",
                    }
                },
            )  # Save as SVG due to issues with PNG saving
            print(f"Visualization saved to: {output_path}")

    except Exception as e:
        logger.error(f"Error saving visualization: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
