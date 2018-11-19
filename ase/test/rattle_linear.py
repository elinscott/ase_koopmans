"""Test RATTLE and QM/MM for rigid linear acetonitrile."""

import numpy as np

from ase import Atoms
from ase.calculators.acn import (ACN, m_me, m_c, 
                                 m_n, r_mec, r_cn,
                                 sigma_me, sigma_c, sigma_n, 
                                 epsilon_me, epsilon_c, epsilon_n)
from ase.calculators.qmmm import SimpleQMMM, EIQMMM, LJInteractionsGeneral
from ase.md.verlet import VelocityVerlet
from ase.constraints import FixBondLengthsLinear 
import ase.units as units

sigma = np.array([sigma_me, sigma_c, sigma_n])
epsilon = np.array([epsilon_me, epsilon_c, epsilon_n])
i = LJInteractionsGeneral(sigma, epsilon, sigma, epsilon)

for calc in [ACN(),
             SimpleQMMM([0, 1, 2], ACN(), ACN(), ACN()),
             EIQMMM([0, 1, 2], ACN(), ACN(), i)]:

    dimer = Atoms('CCNCCN',
                  [(-r_mec, 0, 0),
                   (0, 0, 0),
                   (r_cn, 0, 0),
                   (r_mec, 3.7, 0),
                   (0, 3.7, 0),
                   (-r_cn, 3.7, 0)])

    masses = dimer.get_masses()
    masses[::3] = m_me
    dimer.set_masses(masses) 

    fixd = FixBondLengthsLinear(pairs=[(0, 2), (3, 5)], 
                                singlets=[1, 4],
                                distances=[r_mec,r_cn],
                                masses=[m_me,m_c,m_n]) 
    dimer.set_constraint(fixd)

    dimer.calc = calc

    d1 = dimer[:3].get_all_distances()
    d2 = dimer[3:].get_all_distances()
    e = dimer.get_potential_energy()

    md = VelocityVerlet(dimer, 2.0 * units.fs,
                        trajectory=calc.name + '.traj',
                        logfile=calc.name + '.log',
                        loginterval=20)
    md.run(100)

    de = dimer.get_potential_energy() - e
    
    assert np.all(abs(dimer[:3].get_all_distances()-d1) < 1e-10)
    assert np.all(abs(dimer[3:].get_all_distances()-d2) < 1e-10)
    assert abs(de - -0.005) < 0.001
