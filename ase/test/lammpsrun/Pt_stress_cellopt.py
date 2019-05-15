import numpy as np
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
err = np.abs(stress1 - stress).max()
assert err < 1e-4, 'Stress not as accurate as expected: {}'.format(err)

opt = BFGS(ExpCellFilter(atoms), trajectory='opt.traj')
for i, _ in enumerate(opt.irun(fmax=0.05)):
    pass

assert i < 80, 'Expected 59 iterations, got many more: {}'.format(i)
