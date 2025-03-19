from ase_koopmans.calculators.emt import EMT
from ase_koopmans.io import read
from ase_koopmans.neb import NEB
from ase_koopmans.optimize import BFGS

# read the last structures (of 5 images used in NEB)
images = read('neb.traj@-5:')

for i in range(1, len(images) - 1):
    images[i].calc = EMT()

neb = NEB(images)
qn = BFGS(neb, trajectory='neb_restart.traj')
qn.run(fmax=0.005)
