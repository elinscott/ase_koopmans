def test_strain():
    from math import sqrt
    from ase_koopmans import Atoms
    from ase_koopmans.constraints import StrainFilter
    from ase_koopmans.optimize.mdmin import MDMin
    from ase_koopmans.io import Trajectory
    try:
        from asap3 import EMT
    except ImportError:
        pass
    else:
        a = 3.6
        b = a / 2
        cu = Atoms('Cu', cell=[(0,b,b),(b,0,b),(b,b,0)], pbc=1) * (6, 6, 6)

        cu.calc = EMT()
        f = StrainFilter(cu, [1, 1, 1, 0, 0, 0])
        opt = MDMin(f, dt=0.01)
        t = Trajectory('Cu.traj', 'w', cu)
        opt.attach(t)
        opt.run(0.001)

    # HCP:
        from ase_koopmans.build import bulk
        cu = bulk('Cu', 'hcp', a=a / sqrt(2))
        cu.cell[1,0] -= 0.05
        cu *= (6, 6, 3)

        cu.calc = EMT()
        f = StrainFilter(cu)
        opt = MDMin(f, dt=0.01)
        t = Trajectory('Cu.traj', 'w', cu)
        opt.attach(t)
        opt.run(0.01)

