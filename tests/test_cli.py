"""Tests for the xyz2graph command-line interface."""

import logging
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from xyz2graph import cli
from xyz2graph.graph import MolGraph

from .conftest import WATER_XYZ


@pytest.fixture
def water_file(tmp_path: Path) -> Path:
    path = tmp_path / "water.xyz"
    path.write_text(WATER_XYZ)
    return path


def test_parse_remove_arg_indices_and_elements() -> None:
    indices, elements = cli.parse_remove_arg("1,2,H,O,5")
    assert indices == [1, 2, 5]
    assert elements == ["H", "O"]


def test_parse_remove_arg_skips_blanks_and_invalid(caplog: pytest.LogCaptureFixture) -> None:
    with caplog.at_level(logging.WARNING, logger="xyz2graph"):
        indices, elements = cli.parse_remove_arg("1, ,7,??,O")

    assert indices == [1, 7]
    assert elements == ["O"]
    assert "Ignoring invalid value: ??" in caplog.text


@pytest.mark.parametrize(
    ("indices", "elements", "expected"),
    [
        ([1, 2], [], "Removing atoms at indices: 1, 2"),
        ([], ["H", "O"], "Removing elements: H, O"),
        ([3], ["H"], "Removing atoms at indices: 3; Removing elements: H"),
        ([], [], ""),
    ],
)
def test_format_remove_message(indices: list[int], elements: list[str], expected: str) -> None:
    assert cli.format_remove_message(indices, elements) == expected


def test_generate_output_path_uses_explicit_argument(tmp_path: Path) -> None:
    out = cli.generate_output_path(str(tmp_path / "in.xyz"), str(tmp_path / "out.html"))
    assert out == tmp_path / "out.html"


def test_generate_output_path_falls_back_to_input_stem(tmp_path: Path) -> None:
    out = cli.generate_output_path(str(tmp_path / "in.xyz"), None)
    assert out == tmp_path / "in.html"


def test_log_molecule_state_emits_formula(
    water_molecule: MolGraph, caplog: pytest.LogCaptureFixture
) -> None:
    with caplog.at_level(logging.INFO, logger="xyz2graph"):
        cli.log_molecule_state(water_molecule, "Initial molecule")
    assert "Initial molecule" in caplog.text
    assert "H2O" in caplog.text


def test_log_molecule_state_debug_lists_atoms(
    water_molecule: MolGraph, caplog: pytest.LogCaptureFixture
) -> None:
    cli.logger.setLevel(logging.DEBUG)
    try:
        with caplog.at_level(logging.DEBUG, logger="xyz2graph"):
            cli.log_molecule_state(water_molecule)
        assert "0. O" in caplog.text
    finally:
        cli.logger.setLevel(logging.INFO)


def test_parse_args_minimal(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(sys, "argv", ["xyz2graph", "molecule.xyz"])
    args = cli.parse_args()
    assert args.xyz_file == "molecule.xyz"
    assert args.output is None
    assert args.browser is False
    assert args.debug is False
    assert args.remove is None


def test_parse_args_full(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        sys,
        "argv",
        ["xyz2graph", "mol.xyz", "-o", "out.html", "--browser", "--debug", "-r", "1,H"],
    )
    args = cli.parse_args()
    assert args.output == "out.html"
    assert args.browser is True
    assert args.debug is True
    assert args.remove == "1,H"


def test_main_writes_html_file(
    water_file: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    out = tmp_path / "result.html"
    monkeypatch.setattr(
        sys, "argv", ["xyz2graph", str(water_file), "-o", str(out), "--debug", "-r", "H"]
    )

    cli.main()

    assert out.exists()


def test_main_browser_mode(water_file: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(sys, "argv", ["xyz2graph", str(water_file), "--browser"])
    fake_plot = MagicMock()
    with patch("xyz2graph.cli.offline.plot", fake_plot):
        cli.main()

    assert fake_plot.called


def test_main_handles_missing_file(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(sys, "argv", ["xyz2graph", str(tmp_path / "missing.xyz")])
    with pytest.raises(SystemExit) as exc_info:
        cli.main()
    assert exc_info.value.code == 1


def test_main_handles_save_error(
    water_file: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    out = tmp_path / "out.html"
    monkeypatch.setattr(sys, "argv", ["xyz2graph", str(water_file), "-o", str(out)])

    with patch("xyz2graph.cli.write_html", side_effect=RuntimeError("boom")):
        with pytest.raises(SystemExit) as exc_info:
            cli.main()
    assert exc_info.value.code == 1
