import re
from itertools import combinations
from math import sqrt

atomic_radii = dict(Ac=1.88, Ag=1.59, Al=1.35, Am=1.51, As=1.21, Au=1.50, B=0.83, Ba=1.34, Be=0.35, Bi=1.54, Br=1.21,
                    C=0.68, Ca=0.99, Cd=1.69, Ce=1.83, Cl=0.99, Co=1.33, Cr=1.35, Cs=1.67, Cu=1.52, D=0.23, Dy=1.75,
                    Er=1.73, Eu=1.99, F=0.64, Fe=1.34, Ga=1.22, Gd=1.79, Ge=1.17, H=0.23, Hf=1.57, Hg=1.70, Ho=1.74,
                    I=1.40, In=1.63, Ir=1.32, K=1.33, La=1.87, Li=0.68, Lu=1.72, Mg=1.10, Mn=1.35, Mo=1.47, N=0.68,
                    Na=0.97, Nb=1.48, Nd=1.81, Ni=1.50, Np=1.55, O=0.68, Os=1.37, P=1.05, Pa=1.61, Pb=1.54, Pd=1.50,
                    Pm=1.80, Po=1.68, Pr=1.82, Pt=1.50, Pu=1.53, Ra=1.90, Rb=1.47, Re=1.35, Rh=1.45, Ru=1.40, S=1.02,
                    Sb=1.46, Sc=1.44, Se=1.22, Si=1.20, Sm=1.80, Sn=1.46, Sr=1.12, Ta=1.43, Tb=1.76, Tc=1.35, Te=1.47,
                    Th=1.79, Ti=1.47, Tl=1.55, Tm=1.72, U=1.58, V=1.33, W=1.37, Y=1.78, Yb=1.94, Zn=1.45, Zr=1.56)


class MolGraph:
    """Represents a molecular graph."""
    __slots__ = ['elements', 'x', 'y', 'z', 'adj_list',
                 'atomic_radii', 'bond_lengths']

    def __init__(self):
        self.elements = []
        self.x = []
        self.y = []
        self.z = []
        self.adj_list = {}
        self.atomic_radii = []
        self.bond_lengths = {}

    def read_xyz(self, file_path: str) -> None:
        """Reads an XYZ file, searches for elements and their cartesian coordinates
        and adds them to corresponding arrays."""
        pattern = re.compile(r'([A-Za-z]{1,3})\s*(-?\d+(?:\.\d+)?)\s*(-?\d+(?:\.\d+)?)\s*(-?\d+(?:\.\d+)?)')
        with open(file_path) as file:
            for element, x, y, z in pattern.findall(file.read()):
                self.elements.append(element)
                self.x.append(float(x))
                self.y.append(float(y))
                self.z.append(float(z))
        self.atomic_radii = [atomic_radii[element] for element in self.elements]
        self._generate_adjacency_list()

    def _generate_adjacency_list(self):
        """Generates an adjacency list from atomic cartesian coordinates."""
        node_ids = range(len(self.elements))
        for i, j in combinations(node_ids, 2):
            x_i, y_i, z_i = self.__getitem__(i)[1]
            x_j, y_j, z_j = self.__getitem__(j)[1]
            distance = sqrt((x_i - x_j) ** 2 + (y_i - y_j) ** 2 + (z_i - z_j) ** 2)
            if 0.1 < distance < (self.atomic_radii[i] + self.atomic_radii[j]) * 1.3:
                self.adj_list.setdefault(i, set()).add(j)
                self.adj_list.setdefault(j, set()).add(i)
                self.bond_lengths[frozenset([i, j])] = round(distance, 5)

    def edges(self):
        """Creates an iterator with all graph edges."""
        edges = set()
        for node, neighbours in self.adj_list.items():
            for neighbour in neighbours:
                edge = frozenset([node, neighbour])
                if edge in edges:
                    continue
                edges.add(edge)
                yield node, neighbour

    def __len__(self):
        return len(self.elements)

    def __getitem__(self, position):
        return self.elements[position], (
            self.x[position], self.y[position], self.z[position])
