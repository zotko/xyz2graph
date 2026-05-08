"""Tests covering remaining branches in xyz2graph.graph."""

import logging
from pathlib import Path
from typing import Any
from unittest.mock import patch

import numpy as np
import pytest
from xyz2graph.graph import MolGraph
from xyz2graph.molecule import Atom


def test_indices_property(water_molecule: MolGraph) -> None:
    assert water_molecule.indices == [0, 1, 2]


def test_adj_list(water_molecule: MolGraph) -> None:
    adj = water_molecule.adj_list
    assert adj[0] == {1, 2}
    assert adj[1] == {0}
    assert adj[2] == {0}


def test_from_string_too_few_lines() -> None:
    mol = MolGraph()
    with pytest.raises(ValueError, match="at least 3 lines"):
        mol.from_string("1\ncomment\n")


def test_from_string_atom_count_mismatch_warns(caplog: pytest.LogCaptureFixture) -> None:
    mol = MolGraph()
    content = "5\ncomment\nH 0 0 0\nH 1 0 0\n"
    with caplog.at_level(logging.WARNING, logger="xyz2graph"):
        mol.from_string(content)
    assert "doesn't match" in caplog.text


def test_read_xyz_missing_file(tmp_path: Path) -> None:
    mol = MolGraph()
    with pytest.raises(FileNotFoundError, match="XYZ file not found"):
        mol.read_xyz(tmp_path / "no_such_file.xyz")


def test_write_xyz_handles_oserror(
    water_molecule: MolGraph, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    target = tmp_path / "out.xyz"
    with patch("pathlib.Path.write_text", side_effect=OSError("disk full")):
        water_molecule.write_xyz(target)

    out = capsys.readouterr().out
    assert "Error writing to file" in out


def test_distance_matrix_memory_error_fallback() -> None:
    mol = MolGraph()
    mol.atoms = [
        Atom("H", 0.0, 0.0, 0.0, 0),
        Atom("H", 1.0, 0.0, 0.0, 1),
        Atom("H", 0.0, 1.0, 0.0, 2),
    ]

    real_einsum = np.einsum
    call_count = {"n": 0}

    def fake_einsum(*args: Any, **kwargs: Any) -> Any:  # noqa: ANN401
        call_count["n"] += 1
        if call_count["n"] == 1:
            raise MemoryError("simulated")
        return real_einsum(*args, **kwargs)

    with patch("xyz2graph.graph.np.einsum", side_effect=fake_einsum):
        dist = mol.distance_matrix()

    assert dist.shape == (3, 3)
    assert np.isclose(dist[0, 1], 1.0)
    assert np.isclose(dist[1, 2], np.sqrt(2))


def test_formula_empty_molecule() -> None:
    assert MolGraph().formula() == ""


def test_remove_unused_elements_warning(
    water_molecule: MolGraph, caplog: pytest.LogCaptureFixture
) -> None:
    with caplog.at_level(logging.WARNING, logger="xyz2graph"):
        water_molecule.remove(elements=["Xe"])
    assert "Element(s) not found" in caplog.text
    assert "Xe" in caplog.text


def test_remove_inplace(water_molecule: MolGraph) -> None:
    result = water_molecule.remove(elements=["H"], inplace=True)
    assert result is None
    assert water_molecule.elements == ["O"]
    assert water_molecule.bonds == []


def test_repr_empty() -> None:
    assert repr(MolGraph()) == "MolGraph(empty)"


def test_repr_populated(water_molecule: MolGraph) -> None:
    text = repr(water_molecule)
    assert "MolGraph(" in text
    assert "H2O" in text
    assert "3 atoms" in text
    assert "2 bonds" in text
