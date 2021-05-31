"""
Testing the orthogonal-cell neighbour list used for topology,
in order to have a reference behaviour that we can compare
to with an updated version of it.
"""

import numpy as np
from ase.build import molecule
from pytest import approx

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


def test_distance_matrix():
    at = molecule("C60")
    ase_dij = at.get_all_distances()

    # old implementation
    xyz = at.get_positions().tolist()
    pairs, distances = Topology.distance_matrix(xyz)

    for i, pair in enumerate(pairs):
        assert ase_dij[pair[0], pair[1]] == approx(distances[i])

    # new implementation
    new_pairs, new_distances = Topology.distance_matrix_ase(xyz)
    assert np.all(new_pairs == pairs)
    assert np.all(new_distances == distances)
