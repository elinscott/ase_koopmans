"""
Run some VASP tests to ensure that the VASP calculator works. This
is conditional on the existence of the VASP_COMMAND or VASP_SCRIPT
environment variables

"""

from ase_koopmans.test.vasp import installed2 as installed
from ase_koopmans import Atoms
from ase_koopmans.calculators.vasp import Vasp2 as Vasp
from ase_koopmans.io import read
import numpy as np

assert installed()


def test_main():
    assert installed()

    # simple test calculation of CO molecule
    d = 1.14
    co = Atoms('CO', positions=[(0, 0, 0), (0, 0, d)],
               pbc=True)
    co.center(vacuum=5.)

    calc = Vasp(xc='PBE',
                prec='Low',
                algo='Fast',
                ismear=0,
                sigma=1.,
                istart=0,
                lwave=False,
                lcharg=False,
                ldipol=True)

    co.calc = calc
    energy = co.get_potential_energy()
    forces = co.get_forces()
    dipole_moment = co.get_dipole_moment()

    # check that parsing of vasprun.xml file works
    conf = read('vasprun.xml')
    assert conf.calc.parameters['kpoints_generation']
    assert conf.calc.parameters['sigma'] == 1.0
    assert conf.calc.parameters['ialgo'] == 68
    assert energy - conf.get_potential_energy() == 0.0
    assert np.allclose(conf.get_forces(), forces)
    assert np.allclose(conf.get_dipole_moment(), dipole_moment, atol=1e-6)

    # Cleanup
    calc.clean()
