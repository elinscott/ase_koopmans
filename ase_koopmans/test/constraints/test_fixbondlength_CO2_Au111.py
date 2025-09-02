from math import pi, cos, sin
import pytest

from ase_koopmans import Atoms
from ase_koopmans.calculators.emt import EMT
from ase_koopmans.constraints import FixBondLengths
from ase_koopmans.optimize import BFGS
from ase_koopmans.build import fcc111, add_adsorbate


@pytest.mark.slow
@pytest.mark.parametrize('wrap', [False, True])
def test_CO2_Au111(wrap):
    zpos = cos(134.3 / 2.0 * pi / 180.0) * 1.197
    xpos = sin(134.3 / 2.0 * pi / 180.0) * 1.19
    co2 = Atoms('COO', positions=[(-xpos + 1.2, 0, -zpos),
                                  (-xpos + 1.2, -1.1, -zpos),
                                  (-xpos + 1.2, 1.1, -zpos)])

    slab = fcc111('Au', size=(2, 2, 4), vacuum=2 * 5, orthogonal=True)
    slab.center()
    add_adsorbate(slab, co2, 1.5, 'bridge')
    slab.set_pbc((True, True, False))
    d0 = co2.get_distance(-3, -2)
    d1 = co2.get_distance(-3, -1)

    calc = EMT()
    slab.calc = calc
    if wrap:
        # Remap into the cell so bond is actually wrapped:
        slab.set_scaled_positions(slab.get_scaled_positions() % 1.0)
    constraint = FixBondLengths([[-3, -2], [-3, -1]])
    slab.set_constraint(constraint)
    dyn = BFGS(slab, trajectory='relax_%d.traj' % wrap)
    dyn.run(fmax=0.05)
    assert abs(slab.get_distance(-3, -2, mic=1) - d0) < 1e-9
    assert abs(slab.get_distance(-3, -1, mic=1) - d1) < 1e-9
