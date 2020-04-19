import numpy as np
import pytest

from ase import Atoms
from ase.io import read
from ase.build import bulk


def test_mustem():
    """Check writing and reading a xtl mustem file."""
    # Reproduce the sto xtl file distributed with muSTEM
    atoms = Atoms(['Sr', 'Ti', 'O', 'O', 'O'],
                  scaled_positions=[[0, 0, 0],
                                    [0.5, 0.5, 0.5],
                                    [0.5, 0.5, 0],
                                    [0.5, 0, 0.5],
                                    [0, 0.5, 0.5]],
                  cell=[3.905, 3.905, 3.905],
                  pbc=True)

    filename = 'sto_mustem.xtl'
    STO_DW_dict = {'Sr': 0.78700E-02, 'O': 0.92750E-02, 'Ti': 0.55700E-02}
    STO_DW_dict_Ti_missing = {key: STO_DW_dict[key] for key in ['Sr', 'O']}

    with pytest.raises(TypeError):
        atoms.write(filename)

    with pytest.raises(ValueError):
        atoms.write(filename, keV=300)

    with pytest.raises(TypeError):
        atoms.write(filename,
                    debye_waller_factors=STO_DW_dict)

    atoms.write(filename, keV=300,
                debye_waller_factors=STO_DW_dict)

    atoms2 = read(filename, format='mustem')
    atoms3 = read(filename)

    for _atoms in [atoms2, atoms3]:
        np.testing.assert_allclose(atoms.positions.sum(),
                                   _atoms.positions.sum())
        np.testing.assert_allclose(atoms.cell.sum(), _atoms.cell.sum())

    with pytest.raises(ValueError):
        # Raise an error if there is a missing key.
        atoms.write(filename, keV=300,
                    debye_waller_factors=STO_DW_dict_Ti_missing)

    atoms.write(filename, keV=300,
                debye_waller_factors=STO_DW_dict,
                occupancy={'Sr': 1.0, 'O': 0.5, 'Ti': 0.9})

    with pytest.raises(ValueError):
        # Raise an error if there is a missing key.
        atoms.write(filename, keV=300,
                    debye_waller_factors=STO_DW_dict,
                    occupancy={'O': 0.5, 'Ti': 0.9})

    with pytest.raises(ValueError):
        # Raise an error if the unit cell is not defined.
        atoms4 = Atoms(['Sr', 'Ti', 'O', 'O', 'O'],
                       positions=[[0, 0, 0],
                                  [0.5, 0.5, 0.5],
                                  [0.5, 0.5, 0],
                                  [0.5, 0, 0.5],
                                  [0, 0.5, 0.5]])
        atoms4.write(filename, keV=300,
                     debye_waller_factors=STO_DW_dict)

    # Setting Debye-Waller factor as float.
    Si_atoms = bulk('Si', cubic=True)

    filename = 'Si100.xtl'
    DW = 0.78700E-02
    Si_atoms.write(filename, keV=300, debye_waller_factors=DW)
    Si_atoms2 = read(filename)

    np.testing.assert_allclose(Si_atoms.positions, Si_atoms2.positions)
    np.testing.assert_allclose(Si_atoms.cell, Si_atoms2.cell)
    np.testing.assert_allclose(Si_atoms2.arrays['occupancy'], np.ones(8))
    np.testing.assert_allclose(Si_atoms2.arrays['debye_waller_factors'],
                               np.ones(8) * DW)
