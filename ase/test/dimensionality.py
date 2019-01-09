from ase.geometry import analyze_dimensionality
import ase.build
from ase.lattice.cubic import FaceCenteredCubic


#2D test
atoms = ase.build.mx2(formula='MoS2', kind='2H', a=3.18, thickness=3.19)
atoms.cell[2,2] = 7
atoms.set_pbc((1, 1, 1))
atoms *= 2

intervals = analyze_dimensionality(atoms, method='RDA')
(score, a, b, hr, h, components, cdim) = intervals[0]
assert hr == '2D'

intervals = analyze_dimensionality(atoms, method='TSA')
(score, a, b, hr, h, components, cdim) = intervals[0]
assert hr == '2D'


#3D test
atoms = FaceCenteredCubic(size=(2,2,2), symbol='Cu', pbc=(1,1,1))

intervals = analyze_dimensionality(atoms, method='RDA')
(score, a, b, hr, h, components, cdim) = intervals[0]
assert hr == '3D'

intervals = analyze_dimensionality(atoms, method='TSA')
(score, a, b, hr, h, components, cdim) = intervals[0]
assert hr == '3D'
