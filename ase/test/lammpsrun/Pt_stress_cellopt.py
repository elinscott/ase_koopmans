import numpy as np
from numpy.testing import assert_allclose
from ase.calculators.lammpsrun import LAMMPS
from ase.build import bulk
from ase.test.eam_pot import Pt_u3
from ase.constraints import ExpCellFilter
from ase.optimize import BFGS

# (For now) reuse eam file stuff from other lammps test:
pot_fn = 'Pt_u3.eam'
f = open(pot_fn, 'w')
f.write(Pt_u3)
f.close()
params = {}
params['pair_style'] = 'eam'
params['pair_coeff'] = ['1 1 {}'.format(pot_fn)]
calc = LAMMPS(specorder=['Pt'], files=[pot_fn], **params)

rng = np.random.RandomState(17)

atoms = bulk('Pt') * (2, 2, 2)
atoms.rattle(stdev=0.1)
atoms.cell += 2 * rng.rand(3, 3)
atoms.calc = calc

stress1 = atoms.get_stress()
stress = calc.calculate_numerical_stress(atoms, d=1e-3)
assert_allclose(stress1, stress, atol=1e-4)

opt = BFGS(ExpCellFilter(atoms), trajectory='opt.traj')
for i, _ in enumerate(opt.irun(fmax=0.05)):
    pass

cell1_ref = np.array([
    [0.16298762, 3.89912471, 3.92825365],
    [4.21007577, 0.63362427, 5.04668170],
    [4.42895706, 3.29171414, 0.44623618]])

stress1_ref = np.array([-3.78807088e-05, -2.91140145e-04, -2.72807454e-04,
                        -6.83136671e-05, -8.60692737e-05, -5.91986279e-05])

assert_allclose(np.asarray(atoms.cell), cell1_ref, atol=1e-14)
assert_allclose(atoms.get_stress(), stress1_ref, atol=1e-14)

assert i < 80, 'Expected 59 iterations, got many more: {}'.format(i)
