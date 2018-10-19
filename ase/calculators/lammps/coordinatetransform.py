import decimal as dec
import numpy as np


class Prism(object):

    def __init__(self, cell, pbc=(True, True, True), digits=10):
        """Create a lammps-style triclinic prism object from a cell

        The main purpose of the prism-object is to create suitable
        string representations of prism limits and atom positions
        within the prism.
        When creating the object, the digits parameter (default set to 10)
        specify the precision to use.
        lammps is picky about stuff being within semi-open intervals,
        e.g. for atom positions (when using create_atom in the in-file),
        x must be within [xlo, xhi).
        """
        a, b, c = cell
        an, bn, cn = [np.linalg.norm(v) for v in cell]

        alpha = np.arccos(np.dot(b, c) / (bn * cn))
        beta = np.arccos(np.dot(a, c) / (an * cn))
        gamma = np.arccos(np.dot(a, b) / (an * bn))

        xhi = an
        xyp = np.cos(gamma) * bn
        yhi = np.sin(gamma) * bn
        xzp = np.cos(beta) * cn
        yzp = (bn * cn * np.cos(alpha) - xyp * xzp) / yhi
        zhi = np.sqrt(cn**2 - xzp**2 - yzp**2)

        # Set precision
        self.car_prec = dec.Decimal('10.0') ** \
            int(np.floor(np.log10(max((xhi, yhi, zhi)))) - digits)
        self.dir_prec = dec.Decimal('10.0') ** (-digits)
        self.acc = float(self.car_prec)
        self.eps = np.finfo(xhi).eps

        # For rotating positions from ase to lammps
        Apre = np.array(((xhi, 0, 0),
                         (xyp, yhi, 0),
                         (xzp, yzp, zhi)))
        self.R = np.dot(np.linalg.inv(cell), Apre)

        # Actual lammps cell may be different from what is used to create R
        def fold(vec, pvec, i):
            p = pvec[i]
            x = vec[i] + 0.5 * p
            n = (np.mod(x, p) - x) / p
            return [float(self.f2qdec(a)) for a in (vec + n * pvec)]

        Apre[1, :] = fold(Apre[1, :], Apre[0, :], 0)
        Apre[2, :] = fold(Apre[2, :], Apre[1, :], 1)
        Apre[2, :] = fold(Apre[2, :], Apre[0, :], 0)

        self.A = Apre
        self.Ainv = np.linalg.inv(self.A)

        if self.is_skewed() and \
                (not (pbc[0] and pbc[1] and pbc[2])):
            raise RuntimeError('Skewed lammps cells MUST have '
                               'PBC == True in all directions!')

    def f2qdec(self, f):
        return dec.Decimal(repr(f)).quantize(self.car_prec, dec.ROUND_DOWN)

    def f2qs(self, f):
        return str(self.f2qdec(f))

    def f2s(self, f):
        return str(dec.Decimal(repr(f)).quantize(self.car_prec,
                                                 dec.ROUND_HALF_EVEN))

    def dir2car(self, v):
        """Direct to cartesian coordinates"""
        return np.dot(v, self.A)

    def car2dir(self, v):
        """Cartesian to direct coordinates"""
        return np.dot(v, self.Ainv)

    def fold_to_str(self, v):
        """Fold a position into the lammps cell (semi open)

        Returns tuple of str.
        """
        # Two-stage fold, first into box, then into semi-open interval
        # (within the given precision).
        d = [x % (1 - self.dir_prec) for x in
             map(dec.Decimal,
                 map(repr, np.mod(self.car2dir(v) + self.eps, 1.0)))]
        return tuple([self.f2qs(x) for x in
                      self.dir2car(list(map(float, d)))])

    def get_lammps_prism(self):
        A = self.A
        return A[0, 0], A[1, 1], A[2, 2], A[1, 0], A[2, 0], A[2, 1]

    def get_lammps_prism_str(self):
        """Return a tuple of strings"""
        p = self.get_lammps_prism()
        return tuple([self.f2s(x) for x in p])

    def positions_to_lammps_strs(self, positions):
        """Rotate an ase-cell position to the lammps cell orientation

        Returns tuple of str.
        """
        rot_positions = np.dot(positions, self.R)
        return [tuple([self.f2s(x) for x in position])
                for position in rot_positions]

    def pos_to_lammps_fold_str(self, position):
        """Rotate and fold an ase-cell position into the lammps cell

        Returns tuple of str.
        """
        return self.fold_to_str(np.dot(position, self.R))

    def is_skewed(self):
        acc = self.acc
        prism = self.get_lammps_prism()
        axy, axz, ayz = [np.abs(x) for x in prism[3:]]
        return (axy >= acc) or (axz >= acc) or (ayz >= acc)
