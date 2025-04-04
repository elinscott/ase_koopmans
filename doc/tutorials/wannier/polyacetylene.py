import numpy as np
from gpaw import GPAW

from ase_koopmans import Atoms
from ase_koopmans.dft.kpoints import monkhorst_pack

kpts = monkhorst_pack((13, 1, 1)) + [1e-5, 0, 0]
calc = GPAW(h=.21, xc='PBE', kpts=kpts, nbands=12, txt='poly.txt',
            eigensolver='cg', convergence={'bands': 9})

CC = 1.38
CH = 1.094
a = 2.45
x = a / 2.
y = np.sqrt(CC**2 - x**2)
atoms = Atoms('C2H2', pbc=(True, False, False), cell=(a, 8., 6.),
              calculator=calc, positions=[[0, 0, 0],
                                          [x, y, 0],
                                          [x, y + CH, 0],
                                          [0, -CH, 0]])
atoms.center()
atoms.get_potential_energy()
calc.write('poly.gpw', mode='all')
