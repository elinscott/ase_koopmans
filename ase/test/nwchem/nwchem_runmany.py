import os
import shutil
from ase.build import molecule
from ase.calculators.nwchem import NWChem
from numpy.testing import assert_allclose


def _run_calc(atoms_in, theory, eref, forces=True, **kwargs):
    atoms = atoms_in.copy()
    calc = NWChem(label=theory, theory=theory, **kwargs)
    atoms.set_calculator(calc)
    assert_allclose(atoms.get_potential_energy(), eref, atol=1e-4, rtol=1e-4)
    if forces:
        assert_allclose(atoms.get_forces(),
                        calc.calculate_numerical_forces(atoms),
                        atol=1e-4, rtol=1e-4)
    shutil.rmtree(theory)
    os.remove(theory + '.nwi')
    os.remove(theory + '.nwo')


def main():
    atoms = molecule('H2O')
    # GTO calculations
    _run_calc(atoms, 'dft', -2051.9802410863354, basis='3-21G')
    _run_calc(atoms, 'scf', -2056.7877421222634, basis='3-21G')
    _run_calc(atoms, 'mp2', -2060.1413846247333, basis='3-21G')
    _run_calc(atoms, 'ccsd', -2060.3418911515882, forces=False, basis='3-21G')
    _run_calc(atoms, 'tce', -2060.319141863451, forces=False, basis='3-21G',
              tce=dict(ccd=''))

    atoms.center(vacuum=3)
    atoms.pbc = True
    # Plane wave calculations
    _run_calc(atoms, 'pspw', -465.48899701912626, forces=False,
              nwpw=dict(tolerances='1e-11 1e-11'))
    _run_calc(atoms, 'band', -465.4889948422154, forces=False,
              nwpw=dict(tolerances='1e-11 1e-11'), memory='1024 mb')
    _run_calc(atoms, 'paw', -2065.975544222318, forces=False,
              nwpw=dict(tolerances='1e-11 1e-11'))


main()
