import numpy as np
from pytest import approx

from pygsm.level_of_theories.ase import geom_to_ase


def test_geom_to_ase():
    numbers = [1, 2, 3, 4]
    xyz = np.arange(12).reshape(4, 3)

    atoms = geom_to_ase(numbers, xyz)

    assert atoms.get_chemical_symbols() == ["H", "He", "Li", "Be"]
    assert atoms.get_positions() == approx(xyz)
