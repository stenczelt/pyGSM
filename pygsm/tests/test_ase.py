import numpy as np
from ase.calculators.lj import LennardJones
from pytest import approx, raises

from pygsm.level_of_theories.ase import ASELoT, geom_to_ase, xyz_to_ase
from pygsm.level_of_theories.base_lot import LoTError

xyz_4x4 = [
    ["H", 1.0, 2.0, 3.0],
    ["He", 11.0, 12.0, 13.0],
    ["Li", 21.0, 32.0, 43.0],
    ["Be", 31.0, 32.0, 33.0],
]


def test_geom_to_ase():
    numbers = [1, 2, 3, 4]
    xyz = np.arange(12).reshape(4, 3)

    atoms = geom_to_ase(numbers, xyz)

    assert atoms.get_chemical_symbols() == ["H", "He", "Li", "Be"]
    assert atoms.get_positions() == approx(xyz)


def test_xyz_to_ase():
    atoms = xyz_to_ase(xyz_4x4)
    assert atoms.get_chemical_symbols() == ["H", "He", "Li", "Be"]
    assert atoms.get_positions() == approx(np.array([x[1:] for x in xyz_4x4]))


def test_ase_lot_from_string():
    lot = ASELoT.from_calculator_string(
        calculator_import="ase.calculators.lj.LennardJones",
        calculator_kwargs=dict(epsilon=1.234),
        geom=xyz_4x4,
    )

    assert isinstance(lot.ase_calculator, LennardJones)
    assert lot.ase_calculator.parameters["epsilon"] == approx(1.234)


def test_ase_lot_error():
    with raises(LoTError, match="ASE-calculator's module is not found.*"):
        _ = ASELoT.from_calculator_string(
            calculator_import="ase.calculators.foo.Dummy", geom=xyz_4x4,
        )

    with raises(LoTError, match="ASE-calculator's class.*not found in module .*"):
        _ = ASELoT.from_calculator_string(
            calculator_import="ase.calculators.lj.Dummy", geom=xyz_4x4,
        )
