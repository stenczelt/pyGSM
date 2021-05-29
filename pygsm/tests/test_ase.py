import numpy as np
from pytest import approx

from pygsm.level_of_theories.ase import geom_to_ase, xyz_to_ase


def test_geom_to_ase():
    numbers = [1, 2, 3, 4]
    xyz = np.arange(12).reshape(4, 3)

    atoms = geom_to_ase(numbers, xyz)

    assert atoms.get_chemical_symbols() == ["H", "He", "Li", "Be"]
    assert atoms.get_positions() == approx(xyz)


def test_xyz_to_ase():
    xyz_4x4 = [
        ["H", 1., 2., 3.],
        ["He", 11., 12., 13.],
        ["Li", 21., 32., 43.],
        ["Be", 31., 32., 33.],
    ]

    atoms = xyz_to_ase(xyz_4x4)
    assert atoms.get_chemical_symbols() == ["H", "He", "Li", "Be"]
    assert atoms.get_positions() == approx(np.array([x[1:] for x in xyz_4x4]))
