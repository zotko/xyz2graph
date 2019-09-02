# xyz2graph

xyz2graph - a molecular graph constructor and visualizer.

This script constructs a molecular graph from atomic coordinates.

This tool relies on the most straightforward method for determination of atomic connectivity in molecules,
which is based on interatomic distances and atomic covalent radii. Each interatomic distance is compared to the sum
of covalent radii r<sub>i</sub> and r<sub>j</sub> of both atoms. If the distance between two atoms is within the range d = 1.3(r<sub>i</sub> + r<sub>j</sub>),
that is the sum of covalent radii plus thirty per cent, then an edge connecting the two nodes is added to the
molecular graph.

This script accepts an XYZ ﬁle with chemical elements and their cartesian coordinates as input.

The `Plotly` package is required for visualisation of the molecular graph.
The 'NetworkX' package is required for converting of the molecular graph into the 'NetworkX' graph.

<p align="center">
  <img src="picture.png">
</p>
