import numpy as np
import pytest

from ase.build import bulk
from ase.calculators.test import gradient_test
from ase.calculators.lj import LennardJones
from ase.constraints import UnitCellFilter, ExpCellFilter


@pytest.fixture
def setup_atoms():
    rng = np.random.RandomState(1)
    a0 = bulk('Cu', cubic=True)

    # perturb the atoms
    s = a0.get_scaled_positions()
    s[:, 0] *= 0.995
    a0.set_scaled_positions(s)

    # perturb the cell
    a0.cell += rng.uniform(-1e-1, 1e-2, size=(3, 3))

    a0.set_calculator(LennardJones())
    return a0


@pytest.mark.slow
def test_unitcellfilter(setup_atoms):
    ucf = UnitCellFilter(setup_atoms)
    f, fn = gradient_test(ucf)
    assert abs(f - fn).max() < 1e-6


@pytest.mark.slow
def test_expcellfilter(setup_atoms):
    ecf = ExpCellFilter(setup_atoms)
    # test all derivatives
    f, fn = gradient_test(ecf)
    assert abs(f - fn).max() < 1e-6
