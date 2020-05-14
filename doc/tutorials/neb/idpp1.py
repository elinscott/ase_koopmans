from ase.build import molecule
from ase.neb import NEB
from ase.calculators.emt import EMT
from ase.optimize.fire import FIRE as QuasiNewton

# Optimise molecule.
initial = molecule('C2H6')
initial.calc = EMT()
relax = QuasiNewton(initial)
relax.run(fmax=0.05)

# Create final state.
final = initial.copy()
final.positions[2:5] = initial.positions[[3, 4, 2]]

# Generate blank images.
images = [initial]

for i in range(9):
    images.append(initial.copy())

for image in images:
    image.calc = EMT()
   
images.append(final)

# Run linear interpolation.
neb = NEB(images)
neb.interpolate()

# Run NEB calculation.
qn = QuasiNewton(neb, trajectory='ethane_linear.traj',
                 logfile='ethane_linear.log')
qn.run(fmax=0.05)
