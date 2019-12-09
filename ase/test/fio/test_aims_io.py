from pathlib import Path
import numpy as np
from ase.io.aims import read_aims as read

parent = Path(__file__).parent
format = "aims"

atoms_string = """
lattice_vector 0.0000000000000000 2.7200000000000002 2.7200000000000002
lattice_vector 2.7200000000000002 0.0000000000000000 2.7200000000000002
lattice_vector 2.7200000000000002 2.7200000000000002 0.0000000000000000
atom -0.010000000000000 0.0272000000000000 0.0272000000000000 Si
atom 1.3600000000000001 1.3600000000000001 1.3600000000000001 Si
"""

file = Path("geometry.in")
new_file = Path("geometry.in.tmp")

with open(str(file), "w") as f:
    f.write(atoms_string)

atoms = read(str(file))
file.unlink()


# check cartesian
def test_cartesian():
    """write cartesian coords and check if structure was preserved"""
    atoms.write(str(new_file), format=format)

    new_atoms = read(str(new_file))

    assert np.allclose(atoms.positions, new_atoms.positions)

    new_file.unlink()


# check scaled
def test_scaled():
    """write fractional coords and check if structure was preserved"""
    atoms.write(str(new_file), format=format, scaled=True, wrap=False)

    new_atoms = read(str(new_file))

    assert np.allclose(atoms.positions, new_atoms.positions), (
        atoms.positions,
        new_atoms.positions,
    )

    new_file.unlink()


if __name__ == "__main__":
    test_cartesian()
    test_scaled()
