import numpy as np
from ase.build import bulk, molecule
from ase.io.jsonio import encode, decode


def assert_equal(atoms1, atoms2):
    assert atoms1 == atoms2
    assert set(atoms1.arrays) == set(atoms2.arrays)
    for name in atoms1.arrays:
        assert np.array_equal(atoms1.arrays[name], atoms2.arrays[name]), name


atoms = bulk('Ti')
print('atoms', atoms)
txt = encode(atoms)
print('encoded', txt)

atoms1 = decode(txt)
print('decoded', atoms1)
txt1 = encode(atoms1)
assert txt == txt1
assert_equal(atoms, atoms1)


BeH = molecule('BeH')
assert BeH.has('initial_magmoms')
new_BeH = decode(encode(BeH))
assert_equal(BeH, new_BeH)
assert new_BeH.has('initial_magmoms')
