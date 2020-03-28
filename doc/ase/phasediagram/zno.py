# creates: zno.png
import numpy as np
import matplotlib.pyplot as plt
from ase.phasediagram import Pourbaix, solvated
refs = solvated('Zn')
print(refs)
refs += [('Zn', 0.0), ('ZnO', -3.323), ('ZnO2(aq)', -2.921)]
pb = Pourbaix(refs, Zn=1, O=1)
print(pb.decompose(1.0, 9.0))
U = np.linspace(-2, 2, 200)
pH = np.linspace(-2, 16, 300)
d, names, text = pb.diagram(U, pH, plot=True, show=False)
print(names, text)
assert len(names) == 7
plt.savefig('zno.png')
