"""Tests for the Atom and Bond dataclasses."""

import numpy as np
from xyz2graph.molecule import Atom, Bond


def test_atom_xyz_property() -> None:
    atom = Atom("H", 1.0, 2.0, 3.0, 0)
    assert np.allclose(atom.xyz, np.array([1.0, 2.0, 3.0]))


def test_atom_distance_to() -> None:
    a = Atom("H", 0.0, 0.0, 0.0, 0)
    b = Atom("H", 3.0, 4.0, 0.0, 1)
    assert a.distance_to(b) == 5.0


def test_atom_repr() -> None:
    atom = Atom("O", 0.5, 1.5, 2.5, 7)
    assert repr(atom) == "7. O (0.5, 1.5, 2.5)"


def test_bond_length_post_init() -> None:
    a = Atom("H", 0.0, 0.0, 0.0, 0)
    b = Atom("H", 1.0, 0.0, 0.0, 1)
    bond = Bond(a, b)
    assert bond.length == 1.0


def test_bond_atoms_and_indices() -> None:
    a = Atom("H", 0.0, 0.0, 0.0, 0)
    b = Atom("O", 0.0, 0.0, 1.0, 1)
    bond = Bond(a, b)

    assert bond.atoms == frozenset([a, b])
    assert bond.indices == frozenset([0, 1])


def test_bond_contains() -> None:
    a = Atom("H", 0.0, 0.0, 0.0, 0)
    b = Atom("H", 1.0, 0.0, 0.0, 1)
    c = Atom("O", 5.0, 5.0, 5.0, 2)
    bond = Bond(a, b)

    assert a in bond
    assert b in bond
    assert c not in bond
