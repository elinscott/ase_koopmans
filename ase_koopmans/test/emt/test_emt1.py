def test_emt1():
    from ase_koopmans import Atoms
    from ase_koopmans.calculators.emt import EMT
    from ase_koopmans.constraints import FixBondLength
    from ase_koopmans.io import Trajectory
    from ase_koopmans.optimize import BFGS

    a = 3.6
    b = a / 2
    cu = Atoms('Cu2Ag',
               positions=[(0, 0, 0),
                          (b, b, 0),
                          (a, a, b)],
               calculator=EMT())
    e0 = cu.get_potential_energy()
    print(e0)

    d0 = cu.get_distance(0, 1)
    cu.set_constraint(FixBondLength(0, 1))
    t = Trajectory('cu2ag.traj', 'w', cu)
    qn = BFGS(cu)
    qn.attach(t.write)


    def f():
        print(cu.get_distance(0, 1))


    qn.attach(f)
    qn.run(fmax=0.001)
    assert abs(cu.get_distance(0, 1) - d0) < 1e-14
