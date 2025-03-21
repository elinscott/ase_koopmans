from math import sqrt
from ase_koopmans import Atoms
from ase_koopmans.optimize import LBFGS
from ase_koopmans.constraints import UnitCellFilter
from ase_koopmans.io import Trajectory
from ase_koopmans.optimize.mdmin import MDMin


# XXX This test should have some assertions!  --askhl
def test_unitcellfilter(asap3):
    a = 3.6
    b = a / 2
    cu = Atoms('Cu',
               cell=[(0, b, b), (b, 0, b), (b, b, 0)],
               pbc=1) * (6, 6, 6)
    cu.calc = asap3.EMT()
    f = UnitCellFilter(cu, [1, 1, 1, 0, 0, 0])
    opt = LBFGS(f)
    t = Trajectory('Cu-fcc.traj', 'w', cu)
    opt.attach(t)
    opt.run(5.0)

    # HCP:
    from ase_koopmans.build import bulk
    cu = bulk('Cu', 'hcp', a=a / sqrt(2))
    cu.cell[1,0] -= 0.05
    cu *= (6, 6, 3)
    cu.calc = asap3.EMT()
    print(cu.get_forces())
    print(cu.get_stress())
    f = UnitCellFilter(cu)
    opt = MDMin(f, dt=0.01)
    t = Trajectory('Cu-hcp.traj', 'w', cu)
    opt.attach(t)
    opt.run(0.2)
