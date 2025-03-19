from ase_koopmans.calculators.emt import EMT
from ase_koopmans.cluster import Octahedron
from ase_koopmans.optimize import BFGS

atoms = Octahedron('Ag', 5, cutoff=2)
atoms.calc = EMT()
opt = BFGS(atoms, trajectory='opt.traj')
opt.run(fmax=0.01)
