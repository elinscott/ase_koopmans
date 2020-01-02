import os
from ase.test import require
from ase.test.dftb import Si_Si_skf
from ase.build import bulk
from ase.calculators.dftb import Dftb
from ase.optimize import QuasiNewton
from ase.constraints import ExpCellFilter

require('dftb')

with open('./Si-Si.skf', 'w') as f:
    f.write(Si_Si_skf)

os.environ['DFTB_PREFIX'] = './'

calc = Dftb(label='dftb',
            kpts=(3,3,3),
            Hamiltonian_SCC='Yes')

atoms = bulk('Si')
atoms.set_calculator(calc)

ecf = ExpCellFilter(atoms)
dyn = QuasiNewton(ecf)
dyn.run(fmax=0.01)

e = atoms.get_potential_energy()
assert abs(e - -73.150819) < 1., e
