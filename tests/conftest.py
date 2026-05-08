"""Shared fixtures for the xyz2graph test suite."""

import pytest
from xyz2graph.graph import MolGraph


WATER_XYZ = """3
Water molecule
O 0.0 0.0 0.0
H 0.757 0.586 0.0
H -0.757 0.586 0.0
"""


@pytest.fixture
def water_molecule() -> MolGraph:
    mol = MolGraph()
    mol.from_string(WATER_XYZ)
    return mol
