"""This module defines a 3D geometry utility."""

from typing import NamedTuple, Type


class Point3D(NamedTuple):
    """Represents a point in 3D space."""

    x: float
    y: float
    z: float

    @classmethod
    def midpoint(
        cls: Type["Point3D"], p1: "Point3D", p2: "Point3D", precision: int = 6
    ) -> "Point3D":
        """Calculate the midpoint between two points.

        Args:
            p1: First point
            p2: Second point
            precision: Number of decimal places for rounding (default: 6)

        Returns:
            A new Point3D representing the midpoint
        """
        return cls(
            x=round((p1.x + p2.x) / 2, precision),
            y=round((p1.y + p2.y) / 2, precision),
            z=round((p1.z + p2.z) / 2, precision),
        )

    def __str__(self) -> str:
        """Return a string representation of the point with 6 decimal precision."""
        return f"Point3D(x={self.x:.6f}, y={self.y:.6f}, z={self.z:.6f})"
