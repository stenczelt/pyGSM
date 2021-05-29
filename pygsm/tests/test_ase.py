import numpy as np
from ase.calculators.lj import LennardJones
from pytest import approx, raises

from pygsm.level_of_theories.ase import ASELoT, geom_to_ase, xyz_to_ase
from pygsm.level_of_theories.base_lot import LoTError
from ase.units import Ha, Bohr

xyz_4x4 = [
    ["H", 1.0, 2.0, 3.0],
    ["He", 4.0, 5.0, 6.0],
    ["Li", 7.0, 8.0, 9.0],
    ["Be", 10.0, 11.0, 12.0],
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

    # run not implemented mode
    lot = ASELoT.from_calculator_string(
        calculator_import="ase.calculators.lj.LennardJones", geom=xyz_4x4,
    )

    with raises(
        NotImplementedError,
        match="Run type energgy is not implemented in the ASE calculator interface",
    ):
        # misspelled energy
        lot.run(xyz_4x4, 0, 0, runtype="energgy")


def test_ase_lot_copy():
    lot = ASELoT.from_calculator_string(
        calculator_import="ase.calculators.lj.LennardJones", geom=xyz_4x4,
    )

    copy_lot = ASELoT.copy(lot, dict())

    # we are NOT making a new instance of the calculator
    assert lot.ase_calculator == copy_lot.ase_calculator

    # options matching
    for key in set(lot.options.keys()).union(copy_lot.options.keys()):
        assert lot.options[key] == copy_lot.options[key]


def test_ase_lot_copy_update():
    lot = ASELoT.from_calculator_string(
        calculator_import="ase.calculators.lj.LennardJones", geom=xyz_4x4,
    )

    copy_lot = ASELoT.copy(lot, dict(ID=1))

    # options matching
    for key in lot.options.keys():
        if key == "ID":
            assert copy_lot.options[key] == 1
        else:
            assert lot.options[key] == copy_lot.options[key]


def test_ase_calculation():
    kw = dict(r0=3, rc=10, sigma=3)
    lot = ASELoT.from_calculator_string(
        calculator_import="ase.calculators.lj.LennardJones",
        geom=xyz_4x4,
        calculator_kwargs=kw,
    )

    # ase ref
    calc = LennardJones(**kw)
    atoms = xyz_to_ase(xyz_4x4)
    atoms.calc = calc

    # not doing gradient if not asked for

    # run the calculation
    lot.run(xyz_4x4, 0, 0)

    assert lot._Energies[(0, 0)][0] * Ha == atoms.get_potential_energy()
    assert lot.Gradients[(0, 0)][0] * Ha / Bohr == approx(-atoms.get_forces())


def test_ase_calculation_nograd():
    # not doing gradient if not asked for
    lot = ASELoT.from_calculator_string(
        calculator_import="ase.calculators.lj.LennardJones", geom=xyz_4x4,
    )

    # run the calculation
    lot.run(xyz_4x4, 0, 0, runtype="energy")

    assert (0, 0) in lot._Energies.keys()
    assert (0, 0) not in lot.Gradients.keys()
