import numpy as np
from ase.geometry import wrap_positions


class Prism:
    """The representation of the unit cell in LAMMPS

    The main purpose of the prism-object is to create suitable
    string representations of prism limits and atom positions
    within the prism.

    :param cell: cell in ase coordinate system
    :param pbc: periodic boundaries
    :param tolerance: precision for skewness test
    """

    # !TODO: derive tolerence from cell-dimensions
    def __init__(self, cell, pbc=(True, True, True), tolerance=1.0e-8):
        # Use LQ decomposition to get the lammps cell
        # ase_cell * R = lammps_cell
        qtrans, ltrans = np.linalg.qr(cell.T, mode='complete')
        self.R = qtrans
        self.lammps_cell = ltrans.T
        self.tolerance = tolerance
        self.pbc = pbc

        lc = self.lammps_cell

        # lammps requires positve values on the diagonal of the
        # triangluar matrix -> mirror if necessary
        for i in range(3):
            if lc[i][i] < 0.0:
                mirror_mat = np.eye(3)
                mirror_mat[i][i] = -1.0
                self.lammps_cell = np.dot(mirror_mat, self.lammps_cell.T).T
                self.R = np.dot(self.R, mirror_mat)

        # lammps minimizes the edge length of the parallelepiped
        self.flip = np.array([abs(lc[1][0]/lc[0][0]) > 0.5 and self.pbc[0],
                              abs(lc[2][0]/lc[0][0]) > 0.5 and self.pbc[0],
                              abs(lc[2][1]/lc[1][1]) > 0.5 and self.pbc[1],
                              ])
        if self.flip[0]:
            lc[1][0] -= np.sign(lc[1][0])*lc[0][0]
        if self.flip[1]:
            lc[2][0] -= np.sign(lc[2][0])*lc[0][0]
        if self.flip[2]:
            lc[2][1] -= np.sign(lc[2][1])*lc[1][1]

    def get_lammps_prism(self):
        """Return into lammps coordination system rotated cell

        :returns: lammps cell
        :rtype: np.array

        """
        return self.lammps_cell[(0, 1, 2, 1, 2, 2), (0, 1, 2, 0, 0, 1)]

    # !TODO: detect flip in lammps
    def update_cell(self, xyz, offdiag):
        """Rotate new lammps cell into ase coordinate system

        :param xyz: dimension on the diagonal
        :param offdiag: off-digonal cell elements
        :returns: ase cell
        :rtype: np.arry

        """
        self.lammps_cell = self.to_cell_matrix(xyz, offdiag)
        lc = self.lammps_cell.copy()

        # reverse flip
        if self.flip[0]:
            lc[1][0] -= np.sign(lc[1][0])*lc[0][0]
        if self.flip[1]:
            lc[2][0] -= np.sign(lc[2][0])*lc[0][0]
        if self.flip[2]:
            lc[2][1] -= np.sign(lc[2][1])*lc[1][1]

        return np.dot(lc, self.R.T)

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
        # !TODO: right eps-limit
        # lammps might not like atoms outside the cell
        return wrap_positions(np.dot(vec, self.R), self.lammps_cell,
                              pbc=self.pbc, eps=1e-18)

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
