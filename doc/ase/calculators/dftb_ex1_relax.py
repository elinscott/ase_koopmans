from ase_koopmans.build import molecule
from ase_koopmans.calculators.dftb import Dftb
from ase_koopmans.io import write
from ase_koopmans.optimize import QuasiNewton

atoms = molecule('H2O')
calc = Dftb(atoms=atoms,
            label='h2o',
            Hamiltonian_MaxAngularMomentum_='',
            Hamiltonian_MaxAngularMomentum_O='p',
            Hamiltonian_MaxAngularMomentum_H='s',
            )
atoms.calc = calc

dyn = QuasiNewton(atoms, trajectory='test.traj')
dyn.run(fmax=0.01)
write('final.xyz', atoms)
