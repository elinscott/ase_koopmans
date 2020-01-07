from ase.calculators.dftb import Dftb
from ase.io import write, read
from ase.build import molecule

atoms = molecule('H2O')
calc = Dftb(atoms=atoms,
            label='h2o',
            Driver_='ConjugateGradient',
            Driver_MaxForceComponent=1e-4,
            Driver_MaxSteps=1000,
            Hamiltonian_MaxAngularMomentum_='',
            Hamiltonian_MaxAngularMomentum_O='p',
            Hamiltonian_MaxAngularMomentum_H='s')
atoms.set_calculator(calc)

calc.calculate(atoms)
final = read('geo_end.gen')
write('final.xyz', final)
