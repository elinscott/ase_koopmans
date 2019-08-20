from ase.md.analysis import DiffusionCoefficient
from ase.atoms import Atoms

# Creating a simple trajectory
# Textbook case. The displacement coefficient should be 0.5 A^2 / fs
a = Atoms('He', positions=[(0, 0, 0)])
traj = [a.copy() for i in range(2)]
traj[1].set_positions([(1, 1, 1)])

timestep = 1 #fs
steps_between_images = 1

diffcoeff = DiffusionCoefficient(traj, timestep, steps_between_images, molecule=True)

ans = diffcoeff.calculate(ignore_n_images=0, number_of_segments=1)
ans_orig = 5.0e-06

eps = 1e-10
assert(abs(ans - ans_orig) < eps)
