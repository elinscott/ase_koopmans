#!/usr/bin/env python3

import os
from numpy.testing import assert_allclose

from ase.build import molecule
from ase.calculators.psi4 import Psi4


def test_main():
    atoms = molecule('H2O')
    atoms.rotate(30, 'x')

    calc = Psi4(basis='3-21G')
    atoms.calc = calc

    # Calculate forces ahead of time, compare against finite difference after
    # checking the psi4-calc.dat file
    atoms.get_forces()
    assert_allclose(atoms.get_potential_energy(), -2056.785854116688,
                    rtol=1e-4, atol=1e-4)

    # Test the reader
    calc2 = Psi4()
    calc2.read('psi4-calc')
    assert_allclose(calc2.results['energy'], atoms.get_potential_energy(),
                    rtol=1e-4, atol=1e-4)
    assert_allclose(calc2.results['forces'], atoms.get_forces(),
                    rtol=1e-4, atol=1e-4)

    # Compare analytical vs numerical forces
    assert_allclose(atoms.get_forces(), calc.calculate_numerical_forces(atoms),
                    rtol=1e-4, atol=1e-4)

    os.remove('psi4-calc.dat')
    # Unfortunately, we can't currently remove timer.dat, because Psi4
    # creates the file after this script exits. Not even atexit works.
    # os.remove('timer.dat')
