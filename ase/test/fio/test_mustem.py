import pytest
from ase import Atoms
from ase.io import read


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

    with pytest.raises(TypeError):
        atoms.write(filename)

    with pytest.raises(ValueError):
        atoms.write(filename, keV=300)

    with pytest.raises(TypeError):
        atoms.write(filename,
                    DW={'Sr': 0.78700E-02, 'O': 0.92750E-02, 'Ti': 0.55700E-02})

    atoms.write(filename, keV=300,
                DW={'Sr': 0.78700E-02, 'O': 0.92750E-02, 'Ti': 0.55700E-02})

    atoms2 = read(filename, format='mustem')
    atoms3 = read(filename)

    for _atoms in [atoms2, atoms3]:
        assert atoms.positions.sum() - _atoms.positions.sum() == 0
        assert atoms.cell.sum() - _atoms.cell.sum() == 0

    with pytest.raises(ValueError):
        # Raise an error if there is a missing key.
        atoms.write(filename, keV=300, DW={'Sr': 0.78700E-02, 'O': 0.92750E-02})

    atoms.write(filename, keV=300,
                DW={'Sr': 0.78700E-02, 'O': 0.92750E-02, 'Ti': 0.55700E-02},
                occupancy={'Sr': 1.0, 'O': 0.5, 'Ti': 0.9})

    with pytest.raises(ValueError):
        # Raise an error if there is a missing key.
        atoms.write(filename, keV=300,
                    DW={'Sr': 0.78700E-02, 'O': 0.92750E-02, 'Ti': 0.55700E-02},
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
                     DW={'Sr': 0.78700E-02, 'O': 0.92750E-02, 'Ti': 0.55700E-02})

    # Setting Debye-Waller factor as float.
    Si_atoms = Atoms(symbols='Si'*8,
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

    filename = 'Si100.xtl'
    Si_atoms.write(filename, keV=300, DW=0.78700E-02)
    Si_atoms2 = read(filename)

    assert Si_atoms.positions.sum() - Si_atoms2.positions.sum() == 0
    assert Si_atoms.cell.sum() - Si_atoms2.cell.sum() == 0