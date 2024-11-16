"""This module defines a 3D geometry utility.

Classes:
    Point3D: A class representing a point in 3D space.

Methods:
    Point3D.midpoint(p1: Point3D, p2: Point3D) -> Point3D:
        Calculate the midpoint between two points.
"""

from typing import NamedTuple, Type


class Point3D(NamedTuple):
    """Represents a point in 3D space."""

    x: float
    y: float
    z: float

    @classmethod
    def midpoint(cls: Type["Point3D"], p1: "Point3D", p2: "Point3D") -> "Point3D":
        """Calculate the midpoint between two points."""
        return cls(x=(p1.x + p2.x) / 2, y=(p1.y + p2.y) / 2, z=(p1.z + p2.z) / 2)
