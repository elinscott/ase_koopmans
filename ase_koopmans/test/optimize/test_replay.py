def test_replay():
    from math import sqrt
    from ase_koopmans import Atoms, Atom
    from ase_koopmans.constraints import FixAtoms
    from ase_koopmans.calculators.emt import EMT
    from ase_koopmans.optimize import QuasiNewton
    from ase_koopmans.io import read
    from ase_koopmans.visualize import view

    # Distance between Cu atoms on a (100) surface:
    d = 3.6 / sqrt(2)
    a = Atoms('Cu',
              positions=[(0, 0, 0)],
              cell=(d, d, 1.0),
              pbc=(True, True, False))
    a *= (2, 2, 1)  # 2x2 (100) surface-cell

    # Approximate height of Ag atom on Cu(100) surfece:
    h0 = 2.0
    a += Atom('Ag', (d / 2, d / 2, h0))

    if 0:
        view(a)

    constraint = FixAtoms(range(len(a) - 1))
    a.calc = EMT()
    a.set_constraint(constraint)
    dyn1 = QuasiNewton(a, trajectory='AgCu1.traj', logfile='AgCu1.log')
    dyn1.run(fmax=0.1)

    a = read('AgCu1.traj')
    a.calc = EMT()
    print(a.constraints)
    dyn2 = QuasiNewton(a, trajectory='AgCu2.traj', logfile='AgCu2.log')
    dyn2.replay_trajectory('AgCu1.traj')
    dyn2.run(fmax=0.01)
