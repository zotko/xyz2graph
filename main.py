"""Molecular graph constructor.

This script constructs a molecular graph from atomic coordinates.

This tool relies on the most straightforward method for determination of atomic connectivity in molecules,
which is based on interatomic distances and atomic covalent radii. Each interatomic distance is compared to the sum
of covalent radii r_i and r_j of both atoms.If the distance between two atoms is within the range d = 1.3(r_i + r_j),
that is the sum of covalent radii plus thirty per cent, then an edge connecting the two nodes is added to the
molecular graph.

This script accepts an XYZ ﬁle with chemical elements and their cartesian coordinates as input.

The `plotly` package is required for visualisation of the the molecular graph.
"""

from sys import argv

from molecular_graph import Graph
from plot import plot

script, file = argv

molecule = Graph()

molecule.read_file(file)

plot(adjacency_list=molecule.adjacency_list, elements=molecule.elements, x_coordinates=molecule.x_coordinates,
     y_coordinates=molecule.y_coordinates, z_coordinates=molecule.z_coordinates)
