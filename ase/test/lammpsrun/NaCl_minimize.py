from ase.calculators.lammpsrun import LAMMPS
from ase.spacegroup import crystal
from ase.data import atomic_numbers,  atomic_masses
from ase.optimize import QuasiNewton
from ase.constraints import UnitCellFilter
import numpy as np
from numpy.testing import assert_allclose


a = 6.15
n = 4
nacl = crystal(['Na', 'Cl'], [(0, 0, 0), (0.5, 0.5, 0.5)], spacegroup=225,
               cellpar=[a, a, a, 90, 90, 90]).repeat((n, n, n))

# Buckingham parameters from
# https://physics.stackexchange.com/questions/250018


pair_style = 'buck/coul/long 12.0'
pair_coeff = ['1 1 3796.9 0.2603 124.90']
pair_coeff += ['2 2 1227.2 0.3214 124.90']
pair_coeff += ['1 2 4117.9 0.3048 0.0']
masses = ['1 {}'.format(atomic_masses[atomic_numbers['Na']]),
          '2 {}'.format(atomic_masses[atomic_numbers['Cl']])]

calc = LAMMPS(specorder=['Na', 'Cl'],
              pair_style=pair_style,
              pair_coeff=pair_coeff,
              masses=masses,
              atom_style='charge',
              kspace_style='pppm 1.0e-5',
              keep_tmp_files=True,
              )

for a in nacl:
    if a.symbol == 'Na':
        a.charge = +1.
    else:
        a.charge = -1.

nacl.set_calculator(calc)
stress1_ref = np.array([-4.41412550e-03, -4.41412550e-03, -4.41412550e-03,
                        -2.14623761e-18, -1.49517432e-17, -3.70086806e-20])

assert_allclose(nacl.get_potential_energy(), -1896.216737561538)
assert_allclose(np.asarray(nacl.cell), 24.6 * np.eye(3), atol=1e-14)
assert_allclose(nacl.get_stress(), stress1_ref, atol=1e-14)

E = nacl.get_potential_energy()

ucf = UnitCellFilter(nacl)
dyn = QuasiNewton(ucf, force_consistent=False)
dyn.run(fmax=1.0E-2)

cell2_ref = np.array([
    [2.48422392e+01, 2.01092771e-18, 8.12710402e-16],
    [2.01092771e-18, 2.48422392e+01, 1.22764266e-16],
    [8.12710402e-16, 1.22764266e-16, 2.48422392e+01]])

stress2_ref = np.array([-1.76806144e-04, -1.76806144e-04, -1.76806144e-04,
                        1.27209825e-17, 1.17268975e-17, 4.45395778e-18])

assert_allclose(nacl.get_potential_energy(), -1897.208861729178)
assert_allclose(np.asarray(nacl.cell), cell2_ref, atol=1e-14)
assert_allclose(nacl.get_stress(), stress2_ref, atol=1e-14)
