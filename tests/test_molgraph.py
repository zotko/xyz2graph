import networkx as nx
import numpy as np
import plotly.graph_objs as go
import pytest

from xyz2graph import MolGraph


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
        assert mol_graph.atomic_radii[0] == MolGraph.REFERENCE_RADII["O"]

    def test_to_plotly(self, mol_graph):
        """Test Plotly figure generation."""
        figure = mol_graph.to_plotly()

        # Basic checks (from Copilot)
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
            MolGraph.CPK_COLORS["O"],
            MolGraph.CPK_COLORS["H"],
            MolGraph.CPK_COLORS["H"],
        ]
        assert list(atom_trace.marker.color) == expected_colors

    def test_to_networkx(self, mol_graph):
        """Test NetworkX graph conversion."""
        graph = mol_graph.to_networkx()

        assert graph is not None
        assert len(graph.nodes) == len(mol_graph)
        assert len(graph.edges) > 0

        # Additional specific checks
        assert isinstance(graph, nx.Graph)
        # Check node attributes
        for i, data in graph.nodes(data=True):
            assert "element" in data
            assert "xyz" in data
            assert data["element"] == mol_graph.elements[i]

    def test_generate_adjacency_list(self, mol_graph):
        """Test adjacency list generation."""
        # Basic checks (from Copilot)
        assert mol_graph.adj_list is not None
        assert mol_graph.adj_matrix is not None
        assert len(mol_graph.bond_lengths) > 0

        # Additional specific checks
        assert mol_graph.adj_list[0] == {1, 2}  # O connected to both H
        assert mol_graph.adj_list[1] == {0}  # H1 connected to O
        assert mol_graph.adj_list[2] == {0}  # H2 connected to O

    def test_edges(self, mol_graph):
        """Test edge iteration."""
        edges = list(mol_graph.edges())

        assert len(edges) > 0
        for edge in edges:
            assert len(edge) == 2

        # Additional specific checks
        assert len(edges) == 2  # Water has 2 bonds
        edge_set = {frozenset(edge) for edge in edges}
        assert frozenset([0, 1]) in edge_set  # O-H1 bond
        assert frozenset([0, 2]) in edge_set  # O-H2 bond

    def test_len_and_getitem(self, mol_graph):
        """Test length and indexing operations."""
        # Length check (from Copilot)
        assert len(mol_graph) == len(mol_graph.elements)

        # Getitem check (from Copilot plus additional)
        element, coords = mol_graph[0]
        assert element in MolGraph.REFERENCE_RADII
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
            mol.atomic_radii = [
                mol.REFERENCE_RADII[element] for element in mol.elements
            ]
