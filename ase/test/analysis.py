#test the geometry.analysis module

from ase.geometry.analysis import *
from ase.build import molecule
import numpy as np

mol = molecule('CH3CH2OH')
ana = Analysis(mol)
assert np.shape(ana.adjacency_matrix[0].todense()) == (9, 9)
for imI in range(len(ana.all_bonds)):
    l1 = sum([len(x) for x in ana.all_bonds[imI]])
    l2 = sum([len(x) for x in ana.unique_bonds[imI]])
    assert l1 == l2 * 2

for imi in range(len(ana.all_angles)):
    l1 = sum([len(x) for x in ana.all_angles[imi]])
    l2 = sum([len(x) for x in ana.unique_angles[imi]])
    assert l1 == l2 * 2

for imi in range(len(ana.all_dihedrals)):
    l1 = sum([len(x) for x in ana.all_dihedrals[imi]])
    l2 = sum([len(x) for x in ana.unique_dihedrals[imi]])
    assert l1 == l2 * 2

assert len(ana.get_angles('C','C','H', unique=False)[0]) == len(ana.get_angles('C','C','H', unique=True)[0])*2


