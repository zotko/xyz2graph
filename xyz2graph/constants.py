"""This module defines default constants for atomic radii and CPK colors.

Constants:
    _DEFAULT_RADII (Dict[str, float]): A dictionary mapping element symbols to their default
        atomic radii in angstroms.
    _DEFAULT_CPK_COLORS (Dict[str, str]): A dictionary mapping element symbols to their default
        CPK colors.
"""

from typing import Dict


_DEFAULT_RADII: Dict[str, float] = {
    "Ac": 1.88,
    "Ag": 1.59,
    "Al": 1.35,
    "Am": 1.51,
    "As": 1.21,
    "Au": 1.50,
    "B": 0.83,
    "Ba": 1.34,
    "Be": 0.35,
    "Bi": 1.54,
    "Br": 1.21,
    "C": 0.68,
    "Ca": 0.99,
    "Cd": 1.69,
    "Ce": 1.83,
    "Cl": 0.99,
    "Co": 1.33,
    "Cr": 1.35,
    "Cs": 1.67,
    "Cu": 1.52,
    "D": 0.23,
    "Dy": 1.75,
    "Er": 1.73,
    "Eu": 1.99,
    "F": 0.64,
    "Fe": 1.34,
    "Ga": 1.22,
    "Gd": 1.79,
    "Ge": 1.17,
    "H": 0.23,
    "Hf": 1.57,
    "Hg": 1.70,
    "Ho": 1.74,
    "I": 1.40,
    "In": 1.63,
    "Ir": 1.32,
    "K": 1.33,
    "La": 1.87,
    "Li": 0.68,
    "Lu": 1.72,
    "Mg": 1.10,
    "Mn": 1.35,
    "Mo": 1.47,
    "N": 0.68,
    "Na": 0.97,
    "Nb": 1.48,
    "Nd": 1.81,
    "Ni": 1.50,
    "Np": 1.55,
    "O": 0.68,
    "Os": 1.37,
    "P": 1.05,
    "Pa": 1.61,
    "Pb": 1.54,
    "Pd": 1.50,
    "Pm": 1.80,
    "Po": 1.68,
    "Pr": 1.82,
    "Pt": 1.50,
    "Pu": 1.53,
    "Ra": 1.90,
    "Rb": 1.47,
    "Re": 1.35,
    "Rh": 1.45,
    "Ru": 1.40,
    "S": 1.02,
    "Sb": 1.46,
    "Sc": 1.44,
    "Se": 1.22,
    "Si": 1.20,
    "Sm": 1.80,
    "Sn": 1.46,
    "Sr": 1.12,
    "Ta": 1.43,
    "Tb": 1.76,
    "Tc": 1.35,
    "Te": 1.47,
    "Th": 1.79,
    "Ti": 1.47,
    "Tl": 1.55,
    "Tm": 1.72,
    "U": 1.58,
    "V": 1.33,
    "W": 1.37,
    "Y": 1.78,
    "Yb": 1.94,
    "Zn": 1.45,
    "Zr": 1.56,
}
_DEFAULT_CPK_COLORS: Dict[str, str] = {
    "Ar": "cyan",
    "B": "salmon",
    "Ba": "darkgreen",
    "Be": "darkgreen",
    "Br": "darkred",
    "C": "black",
    "Ca": "darkgreen",
    "Cl": "green",
    "Cs": "violet",
    "F": "green",
    "Fe": "darkorange",
    "Fr": "violet",
    "H": "white",
    "He": "cyan",
    "I": "darkviolet",
    "K": "violet",
    "Kr": "cyan",
    "Li": "violet",
    "Mg": "darkgreen",
    "N": "blue",
    "Na": "violet",
    "Ne": "cyan",
    "O": "red",
    "P": "orange",
    "Ra": "darkgreen",
    "Rb": "violet",
    "S": "yellow",
    "Sr": "darkgreen",
    "Ti": "gray",
    "Xe": "cyan",
}