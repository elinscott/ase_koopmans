from ase.db.core import object_to_bytes, bytes_to_object
import numpy as np

for obj in [1.0,
            {'a': np.zeros((2, 2), np.float32),
             'b': np.zeros((0, 2), int)},
            ['a', 42, True, None, np.nan, np.inf, 1j]]:
    b = object_to_bytes(obj)
    print(b)
    print(bytes_to_object(b))
