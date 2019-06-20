import numpy as np
from ase.geometry import minkowski_reduce

tol = 1E-14
rng = np.random.RandomState(0)

for i in range(40):
    B = rng.uniform(-1, 1, (3, 3))
    R, H = minkowski_reduce(B)
    assert np.allclose(H @ B, R, atol=tol)

    norms = np.linalg.norm(R, axis=1)
    assert (np.argsort(norms) == range(3)).all()

    # Test idempotency
    _, _H = minkowski_reduce(R)
    assert (_H == np.eye(3).astype(np.int)).all()
