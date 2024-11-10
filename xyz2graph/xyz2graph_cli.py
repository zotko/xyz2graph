#!/usr/bin/env python3
import argparse
import logging
import sys
import tempfile
from pathlib import Path
from typing import Union

import plotly.offline as offline
from plotly.io import write_html

from xyz2graph import MolGraph
from xyz2graph.logging import logger


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

    return parser.parse_args()


def generate_output_path(input_path: str, output_path: Union[str, None]) -> Path:
    """Generate the output file path based on input path."""
    if output_path:
        return Path(output_path)
    return Path(input_path).with_suffix(".html")


def main() -> None:
    """Main function for the command-line interface."""
    args = parse_args()

    if args.debug:
        logger.setLevel(logging.DEBUG)
        logger.debug(f"Arguments: {args}")

    # Create MolGraph instance and read XYZ file
    mg = MolGraph()
    try:
        logger.info(f"Reading XYZ file: {args.xyz_file}")
        mg.read_xyz(args.xyz_file)
    except (FileNotFoundError, ValueError) as e:
        logger.error(f"Error reading XYZ file: {e}")
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
            write_html(fig, str(output_path))
            print(f"Visualization saved to: {output_path}")

    except Exception as e:
        logger.error(f"Error saving visualization: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
