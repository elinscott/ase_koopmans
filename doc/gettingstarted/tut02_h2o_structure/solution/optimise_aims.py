import os

from ase_koopmans import Atoms
from ase_koopmans.calculators.aims import Aims
from ase_koopmans.optimize import BFGS

os.environ['ASE_AIMS_COMMAND'] = 'aims.x'
os.environ['AIMS_SPECIES_DIR'] = '/home/alumne/software/FHIaims/species_defaults/light'

atoms = Atoms('HOH',
              positions=[[0, 0, -1], [0, 1, 0], [0, 0, 1]])

calc = Aims(xc='LDA', compute_forces=True)
atoms.calc = calc

opt = BFGS(atoms, trajectory='opt-aims.traj')
opt.run(fmax=0.05)
