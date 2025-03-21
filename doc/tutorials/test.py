# creates: test.txt
from __future__ import print_function

import sys

import neb1
import neb2
import neb3
import numpy as np

# Monkey-patch view() to avoid ASE-GUI windows popping up:
import ase_koopmans.visualize

ase.visualize.view = lambda *args, **kwargs: None

fd = open('test.txt', 'w')

sys.path.append('selfdiffusion')


e1 = np.ptp([i.get_potential_energy() for i in neb1.images])
assert abs(e1 - 0.111) < 0.002

e2 = np.ptp([i.get_potential_energy() for i in neb2.images])
assert abs(e2 - 0.564) < 0.002

e3 = np.ptp([i.get_potential_energy() for i in neb3.images])
assert abs(e3 - 0.239) < 0.002
print(e1, e2, e3, file=fd)
