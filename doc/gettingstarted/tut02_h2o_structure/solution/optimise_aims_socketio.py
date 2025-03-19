import os

from ase_koopmans import Atoms
from ase_koopmans.calculators.aims import Aims
from ase_koopmans.calculators.socketio import SocketIOCalculator
from ase_koopmans.optimize import BFGS

os.environ['ASE_AIMS_COMMAND'] = 'aims.x'
os.environ['AIMS_SPECIES_DIR'] = '/home/alumne/software/FHIaims/species_defaults/light'

atoms = Atoms('HOH',
              positions=[[0, 0, -1], [0, 1, 0], [0, 0, 1]])
opt = BFGS(atoms, trajectory='opt-aims-socketio.traj')

aims = Aims(xc='LDA',
            compute_forces=True,
            use_pimd_wrapper=('UNIX:mysocket', 31415))

with SocketIOCalculator(aims, unixsocket='mysocket') as calc:
    atoms.calc = calc
    opt.run(fmax=0.05)
