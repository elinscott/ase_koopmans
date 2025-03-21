import numpy as np
import pytest

from ase_koopmans.build import bulk
from ase_koopmans.calculators.lj import LennardJones
from ase_koopmans.optimize.precon import PreconLBFGS, Exp
from ase_koopmans.constraints import UnitCellFilter, ExpCellFilter


@pytest.mark.slow
def test_precon():
    cu0 = bulk("Cu") * (2, 2, 2)
    lj = LennardJones(sigma=cu0.get_distance(0,1))

    cu = cu0.copy()
    cu.set_cell(1.2*cu.get_cell())
    cu.calc = lj
    ucf = UnitCellFilter(cu, constant_volume=True)
    opt = PreconLBFGS(ucf, precon=Exp(mu=1.0, mu_c=1.0))
    opt.run(fmax=1e-3)
    assert abs(np.linalg.det(cu.cell)/np.linalg.det(cu0.cell) - 1.2**3) < 1e-3

    # EcpCellFilter allows relaxing to lower tolerance
    cu = cu0.copy()
    cu.set_cell(1.2*cu.get_cell())
    cu.calc = lj
    ecf = ExpCellFilter(cu, constant_volume=True)
    opt = PreconLBFGS(ecf, precon=Exp(mu=1.0, mu_c=1.0))
    opt.run(fmax=1e-3)
    assert abs(np.linalg.det(cu.cell)/np.linalg.det(cu0.cell) - 1.2**3) < 1e-7
