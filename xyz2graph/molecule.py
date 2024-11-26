"""Classes for representing atoms and bonds in a molecular structure.

This module provides fundamental data structures for molecular modeling
"""

from dataclasses import dataclass, field
from typing import FrozenSet, cast

import numpy as np
from numpy.typing import NDArray


@dataclass
class Atom:
    """Represents an atom in a molecular structure."""

    element: str
    x: float
    y: float
    z: float
    index: int
    radius: float = 0.0

    @property
    def xyz(self) -> NDArray[np.float64]:
        """Get atom coordinates as numpy array."""
        return np.array([self.x, self.y, self.z], dtype=np.float64)

    def distance_to(self, other: "Atom") -> float:
        """Calculate distance to another atom."""
        return cast(float, np.linalg.norm(self.xyz - other.xyz).item())

    def __repr__(self) -> str:
        """Return string representation of atom."""
        return f"{self.index}. {self.element} ({self.x}, {self.y}, {self.z})"


@dataclass
class Bond:
    """Represents a chemical bond between two atoms."""

    atom1: Atom
    atom2: Atom
    length: float = field(init=False)

    def __post_init__(self) -> None:
        """Calculate bond length after initialization."""
        dx = self.atom1.x - self.atom2.x
        dy = self.atom1.y - self.atom2.y
        dz = self.atom1.z - self.atom2.z
        self.length = round(np.sqrt(dx * dx + dy * dy + dz * dz), 5)

    @property
    def atoms(self) -> FrozenSet[Atom]:
        """Get the atoms involved in the bond as an order-independent set."""
        return frozenset([self.atom1, self.atom2])

    @property
    def indices(self) -> FrozenSet[int]:
        """Get indices of bonded atoms."""
        return frozenset([self.atom1.index, self.atom2.index])

    def __contains__(self, atom: Atom) -> bool:
        """Check if an atom is part of this bond."""
        return atom in self.atoms
