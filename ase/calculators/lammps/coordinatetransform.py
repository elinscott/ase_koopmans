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

        Qtrans, Ltrans = np.linalg.qr(cell.T)
        self.R = Qtrans
        self.lammps_cell = Ltrans.T
        # Set precision
        self.car_prec = dec.Decimal('10.0') ** \
            int(np.floor(np.log10(np.max(self.lammps_cell)))-digits)
        self.dir_prec = dec.Decimal('10.0') ** (-digits)
        self.acc = float(self.car_prec)
        self.eps = np.finfo(self.acc).eps

        if self.is_skewed() and not np.all(pbc):
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
        return np.dot(v, self.R)

    def car2dir(self, v):
        """Cartesian to direct coordinates"""
        return np.dot(v, self.R.T)

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
        return self.lammps_cell[(0,1,2,1,2,2), (0,1,2,0,0,1)]    

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
        cell_sq = self.lammps_cell**2
        return np.sum(np.tril(cell_sq, -1)) / np.sum(np.diag(cell_sq)) > acc
