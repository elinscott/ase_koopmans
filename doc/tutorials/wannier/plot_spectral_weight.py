import matplotlib.pyplot as plt
import numpy as np
from gpaw import restart

from ase_koopmans.dft import Wannier

atoms, calc = restart('benzene.gpw', txt=None)
wan = Wannier(nwannier=18, calc=calc, fixedstates=15, file='wan18.pickle')

weight_n = np.sum(abs(wan.V_knw[0])**2, 1)
N = len(weight_n)
F = wan.fixedstates_k[0]
plt.figure(1, figsize=(12, 4))
plt.bar(range(1, N + 1), weight_n, width=0.65, bottom=0,
        color='k', edgecolor='k', linewidth=None,
        align='center', orientation='vertical')
plt.plot([F + 0.5, F + 0.5], [0, 1], 'k--')
plt.axis(xmin=0.32, xmax=N + 1.33, ymin=0, ymax=1)
plt.xlabel('Eigenstate')
plt.ylabel('Projection of wannier functions')
plt.savefig('spectral_weight.png')
plt.show()
