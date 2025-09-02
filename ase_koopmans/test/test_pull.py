def test_pull():
    import numpy as np
    from ase_koopmans import Atoms
    from ase_koopmans.calculators.emt import EMT
    from ase_koopmans.io import Trajectory

    Cu = Atoms('Cu', pbc=(1, 0, 0), calculator=EMT())
    traj = Trajectory('Cu.traj', 'w')
    for a in np.linspace(2.0, 4.0, 20):
        Cu.set_cell([a, 1, 1], scale_atoms=True)
        traj.write(Cu)
