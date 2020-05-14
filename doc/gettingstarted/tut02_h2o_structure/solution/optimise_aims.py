import os
from ase import Atoms
from ase.optimize import BFGS
from ase.calculators.aims import Aims

os.environ['ASE_AIMS_COMMAND'] = 'aims.x'
os.environ['AIMS_SPECIES_DIR'] = '/home/alumne/software/FHIaims/species_defaults/light'

atoms = Atoms('HOH',
              positions=[[0, 0, -1], [0, 1, 0], [0, 0, 1]])

calc = Aims(xc='LDA', compute_forces=True)
atoms.calc = calc

opt = BFGS(atoms, trajectory='opt-aims.traj')
opt.run(fmax=0.05)
