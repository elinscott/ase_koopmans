"""Test serialization of ndarrays and other stuff."""

import numpy as np
import io

from ase.io.jsonio import encode, decode, read_json, write_json


assert decode(encode(np.int64(42))) == 42

c = np.array([0.1j])
assert (decode(encode(c)) == c).all()

fd = io.StringIO()

obj1 = {'hello': 'world'}
write_json(fd, obj1)
fd.seek(0)
obj2 = read_json(fd)

print(obj1)
print(obj2)
