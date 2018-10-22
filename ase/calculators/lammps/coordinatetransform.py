import numpy as np


class Prism:
    """The representation of the unit cell in LAMMPS

    The main purpose of the prism-object is to create suitable
    string representations of prism limits and atom positions
    within the prism.

    :param cell: cell in ase coordinate system
    :param pbc: periodic boundaries
    :param tolerance: precision for skewness test
    """

    def __init__(self, cell, pbc=(True, True, True), tolerance=1.0e-6):
        # Use LQ decomposition to get the lammps cell
        # ase_cell * R = lammps_cell
        qtrans, ltrans = np.linalg.qr(cell.T)
        self.R = qtrans
        self.lammps_cell = ltrans.T
        self.tolerance = tolerance

        if self.is_skewed() and not np.all(pbc):
            raise RuntimeError('Skewed lammps cells MUST have '
                               'PBC == True in all directions!')

    def get_lammps_prism(self):
        """Return lammps cell

        :returns: lammps cell
        :rtype: np.array

        """
        return self.lammps_cell[(0, 1, 2, 1, 2, 2), (0, 1, 2, 0, 0, 1)]

    def update_cell(self, xyz, offdiag):
        """Rotate new lammps cell into ase coordinate system

        :param xyz: dimension on the diagonal
        :param offdiag: off-digonal cell elements
        :returns: ase cell
        :rtype: np.arry

        """
        self.lammps_cell = self.to_cell_matrix(xyz, offdiag)
        return np.dot(self.lammps_cell, self.R.T)

    def to_cell_matrix(self, xyz, offdiag):
        """Build lammps triagonal cell from given parameters

        :param xyz: dimensions on the diagonal
        :param offdiag: off-digonal cell elements
        :returns: lammps cell
        :rtype: np.array

        """
        c_xx, c_yy, c_zz = xyz
        c_xy, c_xz, c_yz = offdiag
        return np.array([[c_xx, 0., 0.], [c_xy, c_yy, 0.],
                         [c_xz, c_yz, c_zz]])

    def vector_to_lammps(self, vec):
        """Rotate vector from ase coordinate system to lammps one

        :param vec: to be rotated ase-vector
        :returns: lammps-vector
        :rtype: np.array

        """
        return np.dot(vec, self.R)

    def vector_to_ase(self, vec):
        """Rotate vector from lammps coordinate system to ase one

        :param vec: to be rotated lammps-vector
        :returns: ase-vector
        :rtype: np.array

        """
        return np.dot(vec, self.R.T)

    def is_skewed(self):
        """Test if a lammps cell is not tetragonal

        :returns: bool
        :rtype: bool

        """
        cell_sq = self.lammps_cell**2
        return np.sum(np.tril(cell_sq, -1)) / np.sum(np.diag(cell_sq)) > self.tolerance
