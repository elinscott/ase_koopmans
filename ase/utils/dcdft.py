import numpy as np

from ase.eos import birchmurnaghan


def delta(v1, B1, Bp1, v2, B2, Bp2, symmetric=True, npoints=100):
    if symmetric:
        va = 0.94 * (v1 + v2) / 2
        vb = 1.06 * (v1 + v2) / 2
    else:
        va = 0.94 * v2
        vb = 1.06 * v2
    dv = (vb - va) / npoints
    v = np.linspace(va + dv / 2, vb - dv / 2, npoints)
    e1 = birchmurnaghan(v, 0.0, B1, Bp1, v1)
    e2 = birchmurnaghan(v, 0.0, B2, Bp2, v2)
    return (((e1 - e2)**2).sum() * dv / (vb - va))**0.5
