import numpy as np
import pytest

from ase_koopmans.units import GPa
from ase_koopmans.build import bulk
from ase_koopmans.calculators.test import gradient_test
from ase_koopmans.calculators.lj import LennardJones
from ase_koopmans.constraints import UnitCellFilter, ExpCellFilter
from ase_koopmans.optimize import FIRE, LBFGSLineSearch


@pytest.mark.slow
def test_unitcellfilterpressure():
    a0 = bulk('Cu', cubic=True)

    # perturb the atoms
    s = a0.get_scaled_positions()
    s[:, 0] *= 0.995
    a0.set_scaled_positions(s)

    # perturb the cell
    a0.cell[...] += np.random.uniform(-1e-2, 1e-2,
                                      size=9).reshape((3,3))

    atoms = a0.copy()
    atoms.calc = LennardJones()
    ucf = UnitCellFilter(atoms, scalar_pressure=10.0*GPa)

    # test all derivatives
    f, fn = gradient_test(ucf)
    assert abs(f - fn).max() < 1e-6

    opt = FIRE(ucf)
    opt.run(1e-3)

    # check pressure is within 0.1 GPa of target
    sigma = atoms.get_stress()/GPa
    pressure = -(sigma[0] + sigma[1] + sigma[2])/3.0
    assert abs(pressure - 10.0) < 0.1


    atoms = a0.copy()
    atoms.calc = LennardJones()
    ecf = ExpCellFilter(atoms, scalar_pressure=10.0*GPa)

    # test all deritatives
    f, fn = gradient_test(ecf)
    assert abs(f - fn).max() < 1e-6

    opt = LBFGSLineSearch(ecf)
    opt.run(1e-3)

    # check pressure is within 0.1 GPa of target
    sigma = atoms.get_stress()/GPa
    pressure = -(sigma[0] + sigma[1] + sigma[2])/3.0
    assert abs(pressure - 10.0) < 0.1
