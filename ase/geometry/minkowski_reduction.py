import itertools
import numpy as np


def gauss(B, hu, hv):
    """Calculate a Gauss-reduced lattice basis."""

    u = np.dot(B.T, hu)
    v = np.dot(B.T, hv)

    max_it = 100000    # in practice this is not exceeded
    for it in range(max_it):

        x = int(round(np.dot(u, v) / np.dot(u, u)))
        hu, hv = hv - x * hu, hu

        u = np.dot(B.T, hu)
        v = np.dot(B.T, hv)
        if np.dot(u, u) >= np.dot(v, v):
            return hv, hu

    raise RuntimeError("Gaussian basis not found after %d iterations" % max_it)


def relevant_vectors_2D(u, v):

    cs = np.array([e for e in itertools.product([-1, 0, 1], repeat=2)])
    vs = np.dot(cs, [u, v])
    indices = np.argsort(np.linalg.norm(vs, axis=1))[:7]
    return vs[indices], cs[indices]


def closest_vector(t0, u, v):

    t = t0
    rs, cs = relevant_vectors_2D(u, v)
    a = np.array([0, 0])

    dprev = float("inf")
    max_it = 100000    # in practice this is not exceeded
    for it in range(max_it):

        ds = np.linalg.norm(rs + t, axis=1)
        index = np.argmin(ds)
        if index == 0 or ds[index] >= dprev:
            return a

        dprev = ds[index]
        r = rs[index]
        kopt = int(round(-np.dot(t, r) / np.dot(r, r)))
        a += kopt * cs[index]
        t = t0 + a[0] * u + a[1] * v

    raise RuntimeError("Closest vector not found after %d iterations" % max_it)


def minkowski_reduce(B):

    """Calculate a Minkowski-reduced lattice basis.  The reduced basis
    has the shortest possible vector lengths and has
    norm(a) <= norm(b) <= norm(c).

    Implements the method described in:

    Low-dimensional Lattice Basis Reduction Revisited
    Nguyen, Phong Q. and Stehlé, Damien,
    ACM Trans. Algorithms 5(4) 46:1--46:48, 2009
    https://doi.org/10.1145/1597036.1597050

    Parameters:

    B: array
        The lattice basis to reduce (in row-vector format).

    Returns:

    R: array
        The reduced lattice basis.
    H: array
        The unimodular matrix transformation (R = H @ B).
    """

    H = np.eye(3).astype(np.int)
    norms = np.linalg.norm(B, axis=1)

    max_it = 100000    # in practice this is not exceeded
    for it in range(max_it):

        # Sort vectors by norm
        indices = np.argsort(norms)
        H = H[indices]

        # Gauss-reduce smallest two vectors
        hw = H[2]
        hu, hv = gauss(B, H[0], H[1])
        H[0] = hu
        H[1] = hv

        H = np.array([hu, hv, hw])
        R = H @ B
        u, v, w = R

        X = u / np.linalg.norm(u)
        Y = v - X * np.dot(v, X)
        Y /= np.linalg.norm(Y)

        # Find closest vector to last element of R
        pu, pv, pw = np.dot(R, np.array([X, Y]).T)
        nb = closest_vector(pw, pu, pv)
        hw = np.dot([nb[0], nb[1], 1], H)

        # Update basis
        H = np.array([hu, hv, hw])
        R = H @ B

        norms = np.diag(np.dot(R, R.T))
        if norms[2] >= norms[1] or (nb == 0).all():
            return R, H

    raise RuntimeError("Minkowski basis not found after %d iterations" % max_it)
