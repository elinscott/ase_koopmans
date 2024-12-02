import pytest
from ase_koopmans import Atoms
from ase_koopmans.units import fs
# Rename the class so import pytest won't search for tests on it:
from ase_koopmans.calculators.test import TestPotential as TstPotential
from ase_koopmans.md import VelocityVerlet
from ase_koopmans.io import Trajectory, read
from ase_koopmans.optimize import QuasiNewton
from ase_koopmans.utils import seterr


@pytest.mark.slow
def test_verlet():
    with seterr(all='raise'):
        a = Atoms('4X',
                  masses=[1, 2, 3, 4],
                  positions=[(0, 0, 0),
                             (1, 0, 0),
                             (0, 1, 0),
                             (0.1, 0.2, 0.7)],
                  calculator=TstPotential())
        print(a.get_forces())
        md = VelocityVerlet(a, timestep=0.5 * fs, logfile='-', loginterval=500)
        traj = Trajectory('4N.traj', 'w', a)
        md.attach(traj.write, 100)
        e0 = a.get_total_energy()
        md.run(steps=10000)
        del traj
        assert abs(read('4N.traj').get_total_energy() - e0) < 0.0001

        qn = QuasiNewton(a)
        qn.run(0.001)
        assert abs(a.get_potential_energy() - 1.0) < 0.000002
