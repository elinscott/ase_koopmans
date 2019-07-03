import numpy as np
from ase import Atoms
from ase.calculators.calculator import kpts2sizeandoffsets as k2so

#### Cubic case ####

# Shape from density
kd = 25 / (2 * np.pi)
a = 6.0
size, offsets = map(tuple, k2so(density=kd,
                                atoms=Atoms(cell=(a, a, a), pbc=True)))
assert (size, offsets) == ((5, 5, 5), (0, 0, 0))

# Shape from size option
a = 6.0
size, offsets = map(tuple, k2so(size=(3, 4, 5),))
assert (size, offsets) == ((3, 4, 5), (0, 0, 0))

# Gamma-centering from density
kd = 24 / (2 * np.pi)
size, offsets = map(tuple, k2so(density=kd,
                                gamma=True,
                                atoms=Atoms(cell=(a, a, a), pbc=True)))
assert (size, offsets) == ((4, 4, 4), (0.125, 0.125, 0.125))

# Gamma-centering from size
size, offsets = map(tuple, k2so(size=(3, 4, 5),
                                gamma=True,
                                atoms=Atoms(cell=(a, a, a), pbc=True)))
assert (size, offsets) == ((3, 4, 5), (0., 0.125, 0.))

# Density with irregular shape
cell = [[2, 1, 0], [1, 2, 2], [-1, 0, 2]]
kd = 3
size, offsets = map(tuple, k2so(density=kd,
                                atoms=Atoms(cell=cell, pbc=True)))
assert (size, offsets) == ((29, 22, 26), (0, 0, 0))

# Set even numbers with density
size, offsets = map(tuple, k2so(density=kd,
                                even=True,
                                atoms=Atoms(cell=cell, pbc=True)))
assert (size, offsets) == ((30, 22, 26), (0, 0, 0))

# Set even numbers and Gamma centre with density
size, offsets = map(tuple, k2so(density=kd,
                                even=True,
                                gamma=True,
                                atoms=Atoms(cell=cell, pbc=True)))
assert (size, offsets) == ((30, 22, 26), (1/60, 1/44, 1/52))

# Set odd with density
# Currently failing: kpts2sizeandoffsets does not act as described
# size, offsets = map(tuple, k2so(density=kd,
#                                 even=False,
#                                 atoms=Atoms(cell=cell, pbc=True)))
# assert (size, offsets) == ((29, 23, 27), (0, 0, 0))

#
# set even with size
# Currently failing: rounding not implemented for size option
# size, offsets = map(tuple, k2so(size=(3, 4, 5), even=True))
# assert (size, offsets) == ((4, 4, 6), (0, 0, 0))

# set odd with size
# Currently failing: rounding not implemented for size option
# size, offsets = map(tuple, k2so(size=(3, 4, 5), even=False))
# assert (size, offsets) == ((3, 5, 5), (0, 0, 0))
