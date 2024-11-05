import networkx as nx
import numpy as np
import plotly.graph_objs as go
import pytest

from xyz2graph import DEFAULT_CPK_COLORS, DEFAULT_RADII, MolGraph


class TestMolGraph:
    @pytest.fixture
    def test_xyz_file(self, tmp_path):
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
    def mol_graph(self, test_xyz_file):
        """Create a MolGraph instance with test data."""
        mol = MolGraph()
        mol.read_xyz(test_xyz_file)
        return mol

    def test_read_xyz(self, mol_graph):
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

    def test_to_plotly(self, mol_graph):
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

    def test_len_and_getitem(self, mol_graph):
        """Test length and indexing operations."""
        # Length check
        assert len(mol_graph) == len(mol_graph.elements)

        # Getitem check
        element, coords = mol_graph[0]
        assert element in DEFAULT_RADII
        assert len(coords) == 3
        assert coords == (mol_graph.x[0], mol_graph.y[0], mol_graph.z[0])

    def test_error_handling(self):
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
