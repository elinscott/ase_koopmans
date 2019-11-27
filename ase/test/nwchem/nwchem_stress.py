import os
import shutil
from ase.build import bulk
from ase.calculators.nwchem import NWChem
from numpy.testing import assert_allclose


def main():
    atoms = bulk('C')

    testname = 'stress_test'

    calc = NWChem(theory='pspw',
                  label=testname,
                  nwpw={'lmbfgs': None,
                        'tolerances': '1e-9 1e-9'},
                  )
    atoms.set_calculator(calc)

    assert_allclose(atoms.get_stress(), calc.calculate_numerical_stress(atoms),
                    atol=1e-3, rtol=1e-3)

    shutil.rmtree(testname)
    os.remove(testname + '.nwi')
    os.remove(testname + '.nwo')


main()
