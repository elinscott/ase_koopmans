import os
from ase.build import bulk
from ase.calculators.dftb import Dftb
from ase.optimize import QuasiNewton
from ase.constraints import ExpCellFilter

p = os.path.dirname(__file__)
os.environ['DFTB_PREFIX'] = p if p else './'

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
