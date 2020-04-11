import os
import tempfile
import numpy as np
import pytest

from ase import Atoms
from ase.io import read


tmpDir = tempfile.gettempdir()


def make_Si100_atoms():
    # Reproduce the SI100.XYZ file distributed with prismatic
    atoms = Atoms(symbols='Si'*8,
                  positions=[[0.0000, 0.0000, 0.0000],
                             [2.7150, 2.7150, 0.0000],
                             [1.3575, 4.0725, 1.3575],
                             [4.0725, 1.3575, 1.3575],
                             [2.7150, 0.0000, 2.7150],
                             [0.0000, 2.7150, 2.7150],
                             [1.3575, 1.3575, 4.0725],
                             [4.0725, 4.0725, 4.0725]],
                  cell=[5.43, 5.43, 5.43],
                  pbc=True)

    return atoms


def make_STO_atoms():
    atoms = Atoms(['Sr', 'Ti', 'O', 'O', 'O'],
                  scaled_positions=[[0, 0, 0],
                                    [0.5, 0.5, 0.5],
                                    [0.5, 0.5, 0],
                                    [0.5, 0, 0.5],
                                    [0, 0.5, 0.5]],
                  cell=[3.905, 3.905, 3.905],
                  pbc=True)

    return atoms

def test_write_read_cycle_xyz_prismatic():
    """Check writing and reading a xtl mustem file."""
    # Reproduce the SI100.XYZ file distributed with prismatic
    atoms = make_Si100_atoms()
    atoms.set_array('occupancy', np.ones_like(atoms.numbers))
    atoms.set_array('debye_waller_factor', np.ones_like(atoms.numbers)*0.076)

    filename = os.path.join(tmpDir, 'SI100.XYZ')
    atoms.write(filename=filename, format='prismatic',
                comments='one unit cell of 100 silicon')

    atoms_loaded = read(filename=filename, format='prismatic')

    np.testing.assert_allclose(atoms.positions, atoms_loaded.positions)
    np.testing.assert_allclose(atoms.cell, atoms_loaded.cell)
    np.testing.assert_allclose(atoms.get_array('occupancy'),
                               atoms_loaded.get_array('occupancy'))
    np.testing.assert_allclose(atoms.get_array('debye_waller_factor'),
                               atoms_loaded.get_array('debye_waller_factor'))


def test_write_error():
    """Check missing parameter when writing xyz prismatic file."""
    atoms_Si100 = make_Si100_atoms()
    atoms_STO = make_STO_atoms()
    filename = os.path.join(tmpDir, 'SI100.XYZ')

    with pytest.raises(ValueError):
        # DW not provided
        atoms_Si100.write(filename, format='prismatic')

    # Write file with DW provided as scalar
    atoms_Si100.write(filename, format='prismatic', DW=0.076)

    # Write file with DW provided as dict
    atoms_Si100.write(filename, format='prismatic', DW={'Si':0.076})

    with pytest.raises(ValueError):
        # DW missing keys
        atoms_STO.write(filename, format='prismatic',
                        DW={'Sr': 0.78700E-02, 'O': 0.92750E-02})

    atoms_STO.write(filename, format='prismatic',
                    DW={'Sr': 0.78700E-02, 'O': 0.92750E-02, 'Ti': 0.55700E-02})
