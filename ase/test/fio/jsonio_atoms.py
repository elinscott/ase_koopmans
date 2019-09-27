from ase.build import bulk
from ase.io.jsonio import encode, decode

atoms = bulk('Ti')
print(atoms)

#txt = encode({1:2, 3:4, 'hello': atoms})
txt = encode(atoms)
print(txt)

atoms1 = decode(txt)
print(atoms1)
txt1 = encode(atoms1)

assert txt == txt1
assert atoms == atoms1
