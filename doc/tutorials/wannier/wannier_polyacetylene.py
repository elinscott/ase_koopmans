from __future__ import print_function

import numpy as np
from gpaw import restart

from ase_koopmans.dft import Wannier

atoms, calc = restart('poly.gpw', txt=None)

# Make wannier functions using (one) extra degree of freedom
wan = Wannier(nwannier=6, calc=calc, fixedenergy=1.5)
wan.localize()
wan.save('poly.pickle')
wan.translate_all_to_cell((2, 0, 0))
for i in range(wan.nwannier):
    wan.write_cube(i, 'polyacetylene_%i.cube' % i)

# Print Kohn-Sham bandstructure
ef = calc.get_fermi_level()
f = open('KSbands.txt', 'w')
for k, kpt_c in enumerate(calc.get_ibz_k_points()):
    for eps in calc.get_eigenvalues(kpt=k):
        print(kpt_c[0], eps - ef, file=f)

# Print Wannier bandstructure
f = open('WANbands.txt', 'w')
for k in np.linspace(-.5, .5, 100):
    for eps in np.linalg.eigvalsh(wan.get_hamiltonian_kpoint([k, 0, 0])).real:
        print(k, eps - ef, file=f)
