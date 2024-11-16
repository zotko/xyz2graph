from pathlib import Path
from typing import List, Tuple

import plotly.graph_objs as go
import pytest

from xyz2graph import MolGraph
from xyz2graph.constants import _DEFAULT_CPK_COLORS, _DEFAULT_RADII


class TestMolGraph:
    @pytest.fixture
    def xyz_files(self, tmp_path: Path) -> List[Tuple[str, str, dict]]:
        """Create test XYZ files with different formats."""
        test_cases = [
            # Standard XYZ format
            (
                """3
                Water molecule
                O    0.0    0.0    0.0
                H    0.757  0.586  0.0
                H   -0.757  0.586  0.0
                """,
                "standard.xyz",
                {"xyz_start": 2, "validate": True},
            ),
            # Custom start line format
            (
                """HEADER
                SOME INFO
                ANOTHER LINE
                O  0.0 0.0 0.0
                H  1.0 0.0 0.0
                H -1.0 0.0 0.0""",
                "custom_start.xyz",
                {"xyz_start": 3},
            ),
            # No header format
            (
                """O  0.0 0.0 0.0
                H  1.0 0.0 0.0
                H -1.0 0.0 0.0""",
                "no_header.xyz",
                {"xyz_start": 0},
            ),
            # Empty comment format
            (
                """3

                O    0.0    0.0    0.0
                H    0.757  0.586  0.0
                H   -0.757  0.586  0.0
                """,
                "empty_comment.xyz",
                {"xyz_start": 2},
            ),
        ]

        # Create all test files
        file_paths = []
        for content, filename, params in test_cases:
            file_path = tmp_path / filename
            file_path.write_text(content)
            file_paths.append((str(file_path), content, params))

        return file_paths

    @pytest.mark.parametrize(
        "format_type", ["standard", "custom_start", "no_header", "empty_comment"]
    )
    def test_read_xyz_formats(
        self, xyz_files: List[Tuple[str, str, dict]], format_type: str
    ) -> None:
        """Test reading XYZ files in different formats."""
        # Find the right test case based on format_type
        file_info = next(x for x in xyz_files if format_type in x[0])
        filepath, _, params = file_info

        mol = MolGraph()
        mol.read_xyz(filepath, **params)

        # Common assertions for all formats
        assert len(mol) == 3
        assert mol.elements == ["O", "H", "H"]
        assert len(mol.atomic_radii) == 3
        assert mol.atomic_radii[0] == _DEFAULT_RADII["O"]

        # Format-specific assertions
        if format_type == "standard":
            assert mol.comment == "Water molecule"
        elif format_type == "empty_comment":
            assert mol.comment == ""
        elif format_type == "no_header":
            assert mol.comment == ""

    def test_to_plotly(self, xyz_files: List[Tuple[str, str, dict]]) -> None:
        """Test Plotly figure generation."""
        # Use standard format file for visualization test
        filepath, _, _ = next(x for x in xyz_files if "standard" in x[0])
        mol = MolGraph()
        mol.read_xyz(filepath)

        figure = mol.to_plotly()

        # Basic checks
        assert isinstance(figure, go.Figure)
        assert len(figure.data) == 2  # atom_trace and bond_trace

        # Detailed trace checks
        atom_trace, bond_trace = figure.data
        assert isinstance(atom_trace, go.Scatter3d)
        assert isinstance(bond_trace, go.Scatter3d)
        assert atom_trace.mode == "markers"
        assert bond_trace.mode == "lines"

        # Check colors
        expected_colors = [
            _DEFAULT_CPK_COLORS["O"],
            _DEFAULT_CPK_COLORS["H"],
            _DEFAULT_CPK_COLORS["H"],
        ]
        assert list(atom_trace.marker.color) == expected_colors

    def test_error_handling(self, tmp_path: Path) -> None:
        """Test various error conditions."""
        mol = MolGraph()

        # Test file not found
        with pytest.raises(FileNotFoundError):
            mol.read_xyz("nonexistent.xyz")

        # Test invalid validation
        invalid_content = """2
        Water molecule
        O    0.0    0.0    0.0
        H    0.757  0.586  0.0
        H   -0.757  0.586  0.0"""
        invalid_file = tmp_path / "invalid.xyz"
        invalid_file.write_text(invalid_content)

        with pytest.raises(ValueError, match="Number of coordinates doesn't match atom count"):
            mol.read_xyz(str(invalid_file), validate=True)

        # Test invalid element
        invalid_element_content = """1
            Test
            Xx    0.0    0.0    0.0"""
        invalid_element_file = tmp_path / "invalid_element.xyz"
        invalid_element_file.write_text(invalid_element_content)

        with pytest.raises(ValueError, match="Unknown element found"):
            mol.read_xyz(str(invalid_element_file))

    def test_len_and_getitem(self, xyz_files: List[Tuple[str, str, dict]]) -> None:
        """Test length and indexing operations."""
        filepath, _, _ = next(x for x in xyz_files if "standard" in x[0])
        mol = MolGraph()
        mol.read_xyz(filepath)

        # Length check
        assert len(mol) == len(mol.elements) == 3

        # Getitem check
        element, coords = mol[0]
        assert element == "O"
        assert len(coords) == 3
        assert coords == (mol.x[0], mol.y[0], mol.z[0])

        # Test out of range
        with pytest.raises(IndexError):
            _ = mol[len(mol)]

    def test_formula(self) -> None:
        """Test the molecular formula generation in Hill notation."""
        mg = MolGraph()

        # Test empty molecule
        assert mg.formula() == ""

        # Helper to set up test molecules
        def setup_molecule(elements: List[str]) -> MolGraph:
            mg.elements = elements
            mg.x = [0.0] * len(elements)  # Dummy coordinates
            mg.y = [0.0] * len(elements)
            mg.z = [0.0] * len(elements)
            return mg

        test_cases = [
            (["C", "C", "H", "H", "H", "H", "H", "H"], "C2H6"),  # Ethane
            (["C", "O", "H", "H", "H", "H"], "CH4O"),  # Methanol
            (["H", "O", "H"], "H2O"),  # Water
            (["Na", "Cl"], "ClNa"),  # Sodium chloride
            (["H"], "H"),  # Single atom
            (["C", "C", "C"], "C3"),  # Pure carbon
        ]

        for elements, expected in test_cases:
            mg = setup_molecule(elements)
            assert mg.formula() == expected

    def test_filter(self) -> None:
        """Test the filter functionality."""
        mol = MolGraph()
        mol.elements = ["C", "H", "H", "H", "O", "H"]
        mol.x = [0.0, 1.0, -0.5, -0.5, -1.0, -1.5]
        mol.y = [0.0, 0.0, 0.866, -0.866, 0.0, 0.0]
        mol.z = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        mol.indices = list(range(6))
        mol.atomic_radii = [mol.default_radii[elem] for elem in mol.elements]
        mol._generate_adjacency_list()

        # Test non-inplace filtering
        no_h = mol.filter(elements=["H"])
        assert no_h is not None  # Type guard for mypy
        assert len(no_h) == 2
        assert set(no_h.elements) == {"C", "O"}
        assert no_h.indices == [0, 4]

        # Test filtering by both indices and elements
        result = mol.filter(indices=[1, 2], elements=["O"])
        assert result is not None  # Type guard for mypy
        assert len(result) == 3
        assert "O" not in result.elements
        assert result.indices == [0, 3, 5]

        # Test inplace filtering
        assert mol.filter(elements=["H"], inplace=True) is None
        assert len(mol) == 2
        assert set(mol.elements) == {"C", "O"}

        # Test invalid filtering
        with pytest.raises(ValueError, match="Cannot filter out all atoms"):
            mol.filter(elements=["C", "O"])  # Would remove all atoms
