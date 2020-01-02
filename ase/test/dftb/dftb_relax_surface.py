import os
from ase.test import require
from ase.test.dftb import Si_Si_skf
from ase.build import diamond100
from ase.calculators.dftb import Dftb
from ase.optimize import BFGS
from ase.constraints import FixAtoms

require('dftb')

with open('./Si-Si.skf', 'w') as f:
    f.write(Si_Si_skf)

os.environ['DFTB_PREFIX'] = './'

calc = Dftb(label='dftb',
            kpts=(2, 2, 1),
            Hamiltonian_SCC='Yes',
            Hamiltonian_Filling='Fermi {',
            Hamiltonian_Filling_empty='Temperature [Kelvin] = 500.0',
            )

a = 5.40632280995384
atoms = diamond100('Si', (1, 1, 6), a=a, vacuum=6., orthogonal=True,
                   periodic=True)
atoms.positions[-2:,2] -= 0.2
atoms.set_constraint(FixAtoms(indices=range(4)))
atoms.set_calculator(calc)

dyn = BFGS(atoms, logfile='-')
dyn.run(fmax=0.1)

e = atoms.get_potential_energy()
assert abs(e - -214.036907) < 1., e
