from ase.md.analysis import DiffusionCoefficient
from ase.atoms import Atoms
from ase.io.trajectory import TrajectoryWriter

# Creating a simple traj file
eps = 1e-10
a = Atoms('N3', [(0, 0, 0), (1, 0, 0), (0, 0, 1)])
writer = TrajectoryWriter('N3.traj', mode='w')
writer.write(a)
a.set_positions([(2, 0, 0), (0, 2, 2), (2, 2, 0)])
writer.write(a)
a.set_positions([(4, 0, 0), (0, 4, 4), (4 ,4, 0)])
writer.write(a)

timestep = 1 #fs
steps_between_images = 1

diffcoeff = DiffusionCoefficient('N3.traj', timestep, steps_between_images)
ans = diffcoeff.calculate()

ans_orig = 3.3333333333333335e-05

assert(abs(ans - ans_orig) < eps)

