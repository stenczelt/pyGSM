"""
Testing the orthogonal-cell neighbour list used for topology,
in order to have a reference behaviour that we can compare
to with an updated version of it.
"""

from ase.build import molecule

from pygsm.coordinate_systems.topology import Topology
from pygsm.utilities.elements import ElementData


def test_build_bonds():
    at = molecule("CH4")

    xyz = at.get_positions().tolist()
    elements = [ElementData.from_atomic_number(int(z)) for z in at.get_atomic_numbers()]
    primitive_indices = range(len(at))

    # default
    bonds = Topology.build_bonds(xyz, elements, primitive_indices)

    for i in range(1, 5):
        assert (0, i) in bonds


def test_ase_bonds():
    at = molecule("CH4")

    bonds = Topology.build_bonds_ase(at)

    for i in range(1, 5):
        assert (0, i) in bonds
