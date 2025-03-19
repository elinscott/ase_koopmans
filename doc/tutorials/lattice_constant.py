# creates: lattice_constant.csv

import numpy as np

from ase_koopmans.build import bulk
from ase_koopmans.calculators.emt import EMT
from ase_koopmans.io import Trajectory, read

a0 = 3.52 / np.sqrt(2)
c0 = np.sqrt(8 / 3.0) * a0


traj = Trajectory('Ni.traj', 'w')


eps = 0.01
for a in a0 * np.linspace(1 - eps, 1 + eps, 3):
    for c in c0 * np.linspace(1 - eps, 1 + eps, 3):
        ni = bulk('Ni', 'hcp', a=a, c=c)
        ni.calc = EMT()
        ni.get_potential_energy()
        traj.write(ni)


configs = read('Ni.traj@:')
energies = [config.get_potential_energy() for config in configs]
a = np.array([config.cell[0, 0] for config in configs])
c = np.array([config.cell[2, 2] for config in configs])

functions = np.array([a**0, a, c, a**2, a * c, c**2])
p = np.linalg.lstsq(functions.T, energies, rcond=-1)[0]

p0 = p[0]
p1 = p[1:3]
p2 = np.array([(2 * p[3], p[4]),
               (p[4], 2 * p[5])])
a0, c0 = np.linalg.solve(p2.T, -p1)

fd = open('lattice_constant.csv', 'w')
fd.write('%.3f, %.3f\n' % (a0, c0))
fd.close()
