from pathlib import Path
from typing import List

import numpy as np
import pytest
from xyz2graph.geometry import Point3D
from xyz2graph.graph import Atom, MolGraph


WATER_XYZ = """3
Water molecule
O 0.0 0.0 0.0
H 0.757 0.586 0.0
H -0.757 0.586 0.0
"""


@pytest.fixture
def water_molecule(tmp_path: Path) -> MolGraph:
    """Create a simple water molecule for testing."""
    file_path = tmp_path / "test.xyz"
    file_path.write_text(WATER_XYZ)
    mol = MolGraph()
    mol.read_xyz(file_path)
    return mol


def test_read_xyz(tmp_path: Path) -> None:
    """Test reading a molecule from an XYZ file."""
    file_path = tmp_path / "test.xyz"
    file_path.write_text(WATER_XYZ)
    mol = MolGraph()
    mol.read_xyz(file_path)

    assert len(mol) == 3
    assert mol.comment == "Water molecule"
    assert mol.elements == ["O", "H", "H"]


def test_molecular_geometry(water_molecule: MolGraph) -> None:
    """Test molecular geometry properties."""
    mol = water_molecule
    assert np.allclose(mol.atoms[0].coords, (0.0, 0.0, 0.0))
    assert np.allclose(mol.atoms[1].coords, (0.757, 0.586, 0.0))
    assert np.allclose(mol.atoms[2].coords, (-0.757, 0.586, 0.0))
    assert len(mol.bonds) == 2


def test_bond_properties(water_molecule: MolGraph) -> None:
    """Test bond properties and lengths."""
    mol = water_molecule
    bond_lengths = mol.bond_lengths

    # Check both O-H bonds exist
    assert any(0 in indices and 1 in indices for indices in bond_lengths)
    assert any(0 in indices and 2 in indices for indices in bond_lengths)
    # No H-H bond
    assert not any(1 in indices and 2 in indices for indices in bond_lengths)

    OH_BOND_MIN, OH_BOND_MAX = 0.95, 1.01  # Typical O-H bond length range in Å
    for length in bond_lengths.values():
        assert OH_BOND_MIN < length < OH_BOND_MAX


def test_adjacency_properties(water_molecule: MolGraph) -> None:
    """Test adjacency matrix and list."""
    mol = water_molecule

    # Test adjacency matrix
    adj_matrix = mol.adj_matrix
    assert adj_matrix.shape == (3, 3)
    assert adj_matrix[0, 1] == adj_matrix[1, 0] == 1  # O-H1 bond
    assert adj_matrix[0, 2] == adj_matrix[2, 0] == 1  # O-H2 bond
    assert adj_matrix[1, 2] == adj_matrix[2, 1] == 0  # No H1-H2 bond


def test_distance_matrix(water_molecule: MolGraph) -> None:
    """Test distance matrix calculation."""
    dist_mat = water_molecule.distance_matrix
    assert dist_mat.shape == (3, 3)
    assert np.allclose(dist_mat[0, 0], 0)  # Self distance
    assert np.allclose(dist_mat[0, 1], dist_mat[1, 0])  # Symmetry
    assert np.allclose(dist_mat[1, 2], 2 * 0.757)  # H-H distance


def test_filtering(water_molecule: MolGraph) -> None:
    """Test atom filtering functionality."""
    # Filter hydrogens
    o_only = water_molecule.filter(elements=["H"])
    assert o_only is not None
    assert len(o_only) == 1
    assert o_only.elements == ["O"]


def test_to_networkx(water_molecule: MolGraph) -> None:
    """Test conversion to NetworkX graph."""
    G = water_molecule.to_networkx()
    assert len(G.nodes) == 3
    assert len(G.edges) == 2

    # Check node attributes
    for n in G.nodes:
        assert "element" in G.nodes[n]
        assert "xyz" in G.nodes[n]

    # Check edge attributes
    for _u, _v, data in G.edges(data=True):
        assert "length" in data


def test_error_handling(tmp_path: Path) -> None:
    """Test error handling for various scenarios."""
    mol = MolGraph()

    with pytest.raises(FileNotFoundError):
        mol.read_xyz(tmp_path / "nonexistent.xyz")

    invalid_xyz = """1
Invalid element
X 0.0 0.0 0.0
"""
    invalid_file = tmp_path / "invalid.xyz"
    invalid_file.write_text(invalid_xyz)
    with pytest.raises(ValueError, match="Unknown element"):
        mol.read_xyz(invalid_file)

    malformed_xyz = """1
Malformed coordinates
H 0.0 missing 0.0
"""
    malformed_file = tmp_path / "malformed.xyz"
    malformed_file.write_text(malformed_xyz)
    with pytest.raises(ValueError, match="Invalid coordinate format"):
        mol.read_xyz(malformed_file)


@pytest.mark.parametrize(
    ("elements", "expected_formula"),
    [
        (["C", "C", "H", "H", "H", "H", "H", "H"], "C2H6"),  # Ethane
        (["C", "O", "H", "H", "H", "H"], "CH4O"),  # Methanol
        (["H", "O", "H"], "H2O"),  # Water
        (["Na", "Cl"], "ClNa"),  # Sodium chloride
        (["H"], "H"),  # Single atom
        (["C", "C", "C"], "C3"),  # Pure carbon
        (["Fe", "C", "O", "O", "O", "O", "O"], "CFeO5"),  # Complex case
    ],
)
def test_formula(elements: List[str], expected_formula: str) -> None:
    """Test molecular formula generation."""
    mol = MolGraph()
    for i, el in enumerate(elements):
        pos = Point3D(x=float(i), y=0.0, z=0.0)
        mol.atoms.append(Atom(element=el, position=pos, index=i))
    assert mol.formula() == expected_formula
