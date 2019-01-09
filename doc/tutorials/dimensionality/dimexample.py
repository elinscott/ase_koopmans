import ase.build
from ase.visualize import view
from ase.geometry import analyze_dimensionality


atoms = ase.build.mx2(formula='MoS2', kind='2H', a=3.18, thickness=3.19)
atoms.cell[2,2] = 7
atoms.set_pbc((1, 1, 1))
atoms *= 3

intervals = analyze_dimensionality(atoms, method='RDA')
(score, a, b, hr, h, components, cdim) = intervals[0]
print(sum([e[0] for e in intervals]))
print(hr, h, score, a, b)

atoms.set_tags(components)
view(atoms)
