from pathlib import Path
from typing import List

import plotly.graph_objs as go
import pytest

from xyz2graph import DEFAULT_CPK_COLORS, DEFAULT_RADII, MolGraph


class TestMolGraph:
    @pytest.fixture
    def test_xyz_file(self, tmp_path: Path) -> str:
        """Create a test XYZ file."""
        content = """3
        Water molecule
        O    0.0    0.0    0.0
        H    0.757  0.586  0.0
        H   -0.757  0.586  0.0
        """
        file_path = tmp_path / "test_data.xyz"
        file_path.write_text(content)
        return str(file_path)

    @pytest.fixture
    def mol_graph(self, test_xyz_file: str) -> MolGraph:
        """Create a MolGraph instance with test data."""
        mol = MolGraph()
        mol.read_xyz(test_xyz_file)
        return mol

    def test_read_xyz(self, mol_graph: MolGraph) -> None:
        """Test reading XYZ file and basic structure."""
        assert len(mol_graph) > 0
        assert (
            len(mol_graph.elements)
            == len(mol_graph.x)
            == len(mol_graph.y)
            == len(mol_graph.z)
        )
        assert mol_graph.atomic_radii is not None
        assert mol_graph.adj_list is not None

        # Additional specific checks
        assert mol_graph.elements == ["O", "H", "H"]
        assert len(mol_graph.atomic_radii) == 3
        assert mol_graph.atomic_radii[0] == DEFAULT_RADII["O"]

    def test_to_plotly(self, mol_graph: MolGraph) -> None:
        """Test Plotly figure generation."""
        figure = mol_graph.to_plotly()

        # Basic checks
        assert figure is not None
        assert len(figure.data) == 2  # atom_trace and bond_trace

        # Additional specific checks
        atom_trace, bond_trace = figure.data
        assert isinstance(atom_trace, go.Scatter3d)
        assert isinstance(bond_trace, go.Scatter3d)
        assert atom_trace.mode == "markers"
        assert bond_trace.mode == "lines"

        # Check colors
        expected_colors = [
            DEFAULT_CPK_COLORS["O"],
            DEFAULT_CPK_COLORS["H"],
            DEFAULT_CPK_COLORS["H"],
        ]
        assert list(atom_trace.marker.color) == expected_colors

    # ... rest of the tests remain the same ...

    def test_len_and_getitem(self, mol_graph: MolGraph) -> None:
        """Test length and indexing operations."""
        # Length check
        assert len(mol_graph) == len(mol_graph.elements)

        # Getitem check
        element, coords = mol_graph[0]
        assert element in DEFAULT_RADII
        assert len(coords) == 3
        assert coords == (mol_graph.x[0], mol_graph.y[0], mol_graph.z[0])

    def test_formula(self) -> None:
        """Test the molecular formula generation in Hill notation"""
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

        # Test carbon-containing molecules
        test_cases = [
            # (elements list, expected formula)
            (["C", "C", "H", "H", "H", "H", "H", "H"], "C2H6"),  # Ethane
            (["C", "O", "H", "H", "H", "H"], "CH4O"),  # Methanol
            (["C", "C", "H", "H", "H", "H", "O", "O"], "C2H4O2"),  # Acetic acid
            # Molecules without carbon (alphabetical ordering)
            (["H", "O", "H"], "H2O"),  # Water
            (["N", "H", "H", "H"], "H3N"),  # Ammonia
            (["H", "H", "S", "O", "O", "O", "O"], "H2O4S"),  # Sulfuric acid
            (["Na", "Cl"], "ClNa"),  # Sodium chloride
            # Edge cases
            (["H"], "H"),  # Single hydrogen
            (["O"], "O"),  # Single oxygen
            (["C"], "C"),  # Single carbon
            (["C", "C", "C"], "C3"),  # Pure carbon
            (["H", "H"], "H2"),  # Hydrogen molecule
        ]

        for elements, expected in test_cases:
            mg = setup_molecule(elements)
            assert (
                mg.formula() == expected
            ), f"Failed for {elements}: expected {expected}, got {mg.formula()}"

    def test_xyz_start_custom(self, tmp_path: Path) -> None:
        """Test reading coordinates with custom start line."""
        xyz_content = """HEADER
        SOME INFO
        ANOTHER LINE
        O  0.0 0.0 0.0
        H  1.0 0.0 0.0
        H -1.0 0.0 0.0"""
        xyz_file = tmp_path / "molecule.xyz"
        xyz_file.write_text(xyz_content)

        mol = MolGraph()
        mol.read_xyz(xyz_file, xyz_start=3)

        assert len(mol.elements) == 3
        assert mol.elements == ["O", "H", "H"]
        assert mol.x == [0.0, 1.0, -1.0]

    def test_validate_with_xyz_start_zero(self, tmp_path: Path) -> None:
        """Test validation is disabled when xyz_start=0."""
        xyz_content = """O  0.0 0.0 0.0
        H  1.0 0.0 0.0
        H -1.0 0.0 0.0"""
        xyz_file = tmp_path / "molecule.xyz"
        xyz_file.write_text(xyz_content)

        mol = MolGraph()
        # Should not raise error even though no atom count is present
        mol.read_xyz(xyz_file, xyz_start=0, validate=True)

        assert len(mol.elements) == 3
        assert mol.elements == ["O", "H", "H"]

    def test_validate_true_invalid(self, tmp_path: Path) -> None:
        """Test validation with mismatching atom count."""
        xyz_content = """2
        Water molecule
        O  0.0 0.0 0.0
        H  1.0 0.0 0.0
        H -1.0 0.0 0.0"""
        xyz_file = tmp_path / "molecule.xyz"
        xyz_file.write_text(xyz_content)

        mol = MolGraph()
        with pytest.raises(
            ValueError, match="Number of coordinates doesn't match atom count"
        ):
            mol.read_xyz(xyz_file, validate=True)

    def test_error_handling(self) -> None:
        """Test error conditions."""
        mol = MolGraph()

        # Test file not found
        with pytest.raises(FileNotFoundError):
            mol.read_xyz("nonexistent.xyz")

        # Test invalid element
        mol.elements = ["Xx"]  # Unknown element
        mol.x = [0.0]
        mol.y = [0.0]
        mol.z = [0.0]

        with pytest.raises((KeyError, ValueError)):
            mol.atomic_radii = [mol.default_radii[element] for element in mol.elements]
