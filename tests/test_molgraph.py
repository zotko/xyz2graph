from pathlib import Path
from typing import List

import numpy as np
import pytest
from xyz2graph.graph import MolGraph
from xyz2graph.molecule import Atom


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
    assert np.allclose(mol.atoms[0].xyz, (0.0, 0.0, 0.0))
    assert np.allclose(mol.atoms[1].xyz, (0.757, 0.586, 0.0))
    assert np.allclose(mol.atoms[2].xyz, (-0.757, 0.586, 0.0))
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
    dist_mat = water_molecule.distance_matrix()
    assert dist_mat.shape == (3, 3)
    assert np.allclose(dist_mat[0, 0], 0)  # Self distance
    assert np.allclose(dist_mat[0, 1], dist_mat[1, 0])  # Symmetry
    assert np.allclose(dist_mat[1, 2], 2 * 0.757)  # H-H distance


def test_remove(water_molecule: MolGraph) -> None:
    """Test atom filtering functionality."""
    # Filter hydrogens
    o_only = water_molecule.remove(elements=["H"])
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


def test_set_element_radius() -> None:
    """Test custom radius setting for elements."""
    mol = MolGraph()
    mol.atoms = [Atom("H", 0, 0, 0, 0), Atom("H", 1, 0, 0, 1)]

    mol.set_element_radius("H", 2.0)
    assert all(atom.radius == 2.0 for atom in mol.atoms)
    assert len(mol.bonds) > 0


def test_invalid_xyz_content(tmp_path: Path) -> None:
    """Test error handling for invalid XYZ file content."""
    invalid_cases = [
        "not_a_number\ncomment\nH 0 0 0",
        "1\ncomment\nH invalid coords",
        "1\ncomment\nXx 0 0 0",  # Invalid element
    ]

    mol = MolGraph()
    for content in invalid_cases:
        xyz_file = tmp_path / "test.xyz"
        xyz_file.write_text(content)
        with pytest.raises((ValueError, FileNotFoundError)):
            mol.read_xyz(xyz_file)


def test_atom_indexing() -> None:
    """Test atom indexing and iteration."""
    mol = MolGraph()
    atoms = [Atom("H", 0, 0, 0, 0), Atom("O", 1, 1, 1, 1)]
    mol.atoms = atoms

    assert mol[0] == atoms[0]
    assert mol[1] == atoms[1]
    with pytest.raises(IndexError):
        _ = mol[2]


def test_remove_invalid_cases() -> None:
    """Test error cases for remove method."""
    mol = MolGraph()
    mol.atoms = [Atom("H", i, 0, 0, i) for i in range(3)]

    with pytest.raises(IndexError, match="Atom index out of range"):
        mol.remove(indices=[10])

    with pytest.raises(ValueError, match="Cannot remove all atoms from molecule"):
        mol.remove(indices=[0, 1, 2])


def test_property_consistency() -> None:
    """Test consistency between different property representations."""
    mol = MolGraph()
    mol.atoms = [Atom("H", 0, 1, 2, 0), Atom("O", 3, 4, 5, 1)]

    assert np.allclose(mol.xyz, np.array([[0, 1, 2], [3, 4, 5]]))
    assert mol.x == [0, 3]
    assert mol.y == [1, 4]
    assert mol.z == [2, 5]


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
        mol.atoms.append(Atom(element=el, x=float(i), y=0.0, z=0.0, index=i))
    assert mol.formula() == expected_formula
