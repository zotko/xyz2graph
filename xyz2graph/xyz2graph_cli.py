#!/usr/bin/env python3
import argparse
import tempfile
from pathlib import Path

import plotly.offline as offline
from plotly.io import write_html

from xyz2graph import MolGraph


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

    return parser.parse_args()


def generate_output_path(input_path: str, output_path: str | None) -> Path:
    """Generate the output file path based on input path."""
    if output_path:
        return Path(output_path)

    input_path = Path(input_path)
    return input_path.with_suffix(".html")


def main() -> None:
    """Main function for the command-line interface."""
    args = parse_args()

    # Create MolGraph instance and read XYZ file
    mg = MolGraph()
    try:
        mg.read_xyz(args.xyz_file)
    except (FileNotFoundError, ValueError) as e:
        print(f"Error: {e}")
        return

    # Generate figure
    fig = mg.to_plotly()

    # Save/display the visualization
    try:
        if args.browser:
            # Create a temporary file that will be automatically cleaned up
            with tempfile.NamedTemporaryFile(
                prefix=Path(args.xyz_file).stem + "_", suffix=".html", delete=False
            ) as tmp:
                offline.plot(fig, filename=tmp.name, auto_open=True)
        else:
            output_path = generate_output_path(args.xyz_file, args.output)
            write_html(fig, str(output_path))
            print(f"Visualization saved to: {output_path}")

    except Exception as e:
        print(f"Error saving visualization: {e}")
        return


if __name__ == "__main__":
    main()
