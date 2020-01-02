import os
from ase import Atoms
from ase.test import require
from ase.test.dftb import Si_Si_skf
from ase.calculators.dftb import Dftb
from ase.optimize import BFGS

require('dftb')

with open('./Si-Si.skf', 'w') as f:
    f.write(Si_Si_skf)

os.environ['DFTB_PREFIX'] = './'

calc = Dftb(label='dftb',
            Hamiltonian_SCC='No',
            Hamiltonian_PolynomialRepulsive='SetForAll {Yes}',
            )

atoms = Atoms('Si2', positions=[[5., 5., 5.], [7., 5., 5.]],
              cell=[12.]*3, pbc=False)
atoms.set_calculator(calc)

dyn = BFGS(atoms, logfile='-')
dyn.run(fmax=0.1)

e = atoms.get_potential_energy()
assert abs(e - -64.830901) < 1., e
