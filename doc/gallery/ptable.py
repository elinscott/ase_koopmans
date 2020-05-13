# creates: ptable.png
import numpy as np
from ase import Atoms

# generates column, row positions for each element
zmax = 118
x, y = 1, 1  # column, row , initial coordinates for Hydrogen
positions = np.zeros((zmax + 1, 3))
for z in range(1, zmax + 1):
    if z == 2:
        x += 16
    if z == 5:
        x += 10
    if z == 13:
        x += 10
    if z == 57 or z == 89:
        y += 3
    if z == 72 or z == 104:
        y -= 3
        x -= 14
    positions[z] = (x, -y, 0)
    x += 1
    if x > 18:
        x = 1
        y += 1

spacing = 2.5
atoms = Atoms(np.arange(1, zmax + 1),
              positions[1:] * spacing)

atoms.write('ptable.png')

# from ase.visualize import view
# view(atoms)
