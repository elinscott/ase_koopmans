import ase.build
from ase.lattice.cubic import FaceCenteredCubic
from ase.geometry.dimensionality import analyze_dimensionality


# 2D test
atoms = ase.build.mx2(formula='MoS2', kind='2H', a=3.18, thickness=3.19)
atoms.cell[2, 2] = 7
atoms.set_pbc((1, 1, 1))
atoms *= 2

intervals = analyze_dimensionality(atoms, method='RDA')
m = intervals[0]
assert m.dimtype == '2D'

intervals = analyze_dimensionality(atoms, method='TSA')
m = intervals[0]
assert m.dimtype == '2D'


# 3D test
atoms = FaceCenteredCubic(size=(2, 2, 2), symbol='Cu', pbc=(1, 1, 1))

intervals = analyze_dimensionality(atoms, method='RDA')
m = intervals[0]
assert m.dimtype == '3D'

intervals = analyze_dimensionality(atoms, method='TSA')
m = intervals[0]
assert m.dimtype == '3D'
