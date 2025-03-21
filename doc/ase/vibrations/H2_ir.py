from gpaw import GPAW, FermiDirac
from gpaw.cluster import Cluster

from ase_koopmans import optimize
from ase_koopmans.build import molecule
from ase_koopmans.vibrations.infrared import InfraRed

h = 0.22

atoms = Cluster(molecule('H2'))
atoms.minimal_box(3.5, h=h)

# relax the molecule
calc = GPAW(h=h, occupations=FermiDirac(width=0.1))
atoms.calc = calc
dyn = optimize.FIRE(atoms)
dyn.run(fmax=0.05)
atoms.write('relaxed.traj')

# finite displacement for vibrations
ir = InfraRed(atoms)
ir.run()
