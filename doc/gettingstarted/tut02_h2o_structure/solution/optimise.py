from gpaw import GPAW

from ase_koopmans import Atoms
from ase_koopmans.optimize import BFGS

atoms = Atoms('HOH',
              positions=[[0, 0, -1], [0, 1, 0], [0, 0, 1]])
atoms.center(vacuum=3.0)

calc = GPAW(mode='lcao', basis='dzp', txt='gpaw.txt')
atoms.calc = calc

opt = BFGS(atoms, trajectory='opt.traj')
opt.run(fmax=0.05)
