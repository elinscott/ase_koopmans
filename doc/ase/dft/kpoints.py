# creates: cc.png
import matplotlib.pyplot as plt
import numpy as np

from ase_koopmans.dft.kpoints import cc162_1x1

B = [(1, 0, 0), (-0.5, 3**0.5 / 2, 0), (0, 0, 1)]
k = np.dot(cc162_1x1, B)
plt.figure(figsize=(5, 4))
plt.plot(k[:, 0], k[:, 1], 'o')
plt.savefig('cc.png')
