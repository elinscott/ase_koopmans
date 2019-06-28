import numpy as np
from ase.utils.arraywrapper import arraylike


__all__ = ['Cell']


@arraylike
class Cell:
    """Parallel epipedal unit cell of up to three dimensions.

    This object resembles a 3x3 array whose [i, j]-th element is the jth
    Cartesian coordinate of the ith unit vector.

    Cells of less than three dimensions are represented by placeholder
    unit vectors that are zero."""

    ase_objtype = 'cell'  # For JSON'ing

    def __init__(self, array, pbc=None):
        if pbc is None:
            # pbc defaults to whether each cell vector is nonzero:
            pbc = array.any(1)

        # We could have lazy attributes for structure (bcc, fcc, ...)
        # and other things.  However this requires making the cell
        # array readonly, else people will modify it and things will
        # be out of synch.
        assert array.shape == (3, 3)
        assert array.dtype == float
        assert pbc.shape == (3,)
        assert pbc.dtype == bool
        self.array = array
        self.pbc = pbc

    def cellpar(self, radians=False):
        """Get cell lengths and angles of this cell.

        See also :func:`ase.geometry.cell.cell_to_cellpar`."""
        from ase.geometry.cell import cell_to_cellpar
        return cell_to_cellpar(self.array, radians)

    def todict(self):
        return dict(array=self.array, pbc=self.pbc)

    @property
    def shape(self):
        return self.array.shape

    @classmethod
    def ascell(cls, cell):
        """Return argument as a Cell object.

        A new Cell object is created if necessary."""
        if isinstance(cell, cls):
            return cell
        return cls.new(cell)

    @classmethod
    def new(cls, cell=None, pbc=None):
        """Create new cell from any parameters."""
        if pbc is None:
            pbc = getattr(cell, 'pbc', None)

        if cell is None:
            cell = np.zeros((3, 3))

        cell = np.array(cell, float)

        if cell.shape == (3,):
            cell = np.diag(cell)
        elif cell.shape == (6,):
            from ase.geometry.cell import cellpar_to_cell
            cell = cellpar_to_cell(cell)
        elif cell.shape != (3, 3):
            raise ValueError('Cell must be length 3 sequence, length 6 '
                             'sequence or 3x3 matrix!')

        cellobj = cls(cell)
        if pbc is not None:
            cellobj.pbc[:] = pbc

        return cellobj

    @classmethod
    def fromcellpar(cls, cellpar, ab_normal=(0, 0, 1), a_direction=None,
                    pbc=None):
        """Return new Cell from cell lengths and angles.

        See also :func:`~ase.geometry.cell.cellpar_to_cell()`."""
        from ase.geometry.cell import cellpar_to_cell
        cell = cellpar_to_cell(cellpar, ab_normal, a_direction)
        return cls(cell, pbc=pbc)

    def get_bravais_lattice(self, eps=2e-4):
        """Return :class:`~ase.lattice.BravaisLattice` for this cell::

        >>> cell = Cell.fromcellpar([4, 4, 4, 60, 60, 60])
        >>> print(cell.get_bravais_lattice())
        FCC(a=5.65685)

        """
        from ase.lattice import get_bravais_lattice
        return get_bravais_lattice(self, eps=eps)

    def bandpath(self, path=None, npoints=None, density=None,
                 special_points=None, eps=2e-4):
        """Build a :class:`~ase.dft.kpoints.BandPath` for this cell.

        If special points are None, determine the Bravais lattice of
        this cell and return a suitable Brillouin zone path with
        standard special points.  With special points, interpolate the
        path directly without determining Bravais lattice.

        Parameters:

        path: string
            String of special point names defining the path, e.g. 'GXL'.
        npoints: int
            Number of points in total.  Note that at least one point
            is added for each special point in the path.
        density: float
            density of kpoints along the path in Å⁻¹.
        special_points: dict
            Dictionary mapping special points to scaled kpoint coordinates.
        eps: float
            Tolerance for determining Bravais lattice.

        Example
        -------
        >>> cell = Cell.fromcellpar([4, 4, 4, 60, 60, 60])
        >>> cell.bandpath('GXW', npoints=20)
        BandPath(path='GXW', special_points={GKLUWX}, kpts=[20 kpoints])

        """
        # TODO: Combine with the rotation transformation from bandpath()
        if special_points is None:
            from ase.lattice import identify_lattice
            lat, op = identify_lattice(self, eps=eps)
            path = lat.bandpath(path, npoints=npoints, density=density)
            return path.transform(op)
        else:
            from ase.dft.kpoints import BandPath
            path = BandPath(path, special_points=special_points)
            return path.interpolate(npoints=npoints, density=density)


    # XXX adapt the transformation stuff and include in the bandpath method.
    def oldbandpath(self, path=None, npoints=None, density=None, eps=2e-4):
        bravais = self.get_bravais_lattice(eps=eps)
        transformation = bravais.get_transformation(self.array)
        return bravais.bandpath(path=path, npoints=npoints, density=density,
                                transformation=transformation)

    def complete(self):
        """Convert missing cell vectors into orthogonal unit vectors."""
        from ase.geometry.cell import complete_cell
        return Cell(complete_cell(self.array), self.pbc.copy())

    def copy(self):
        """Return a copy of this cell."""
        return Cell(self.array.copy(), self.pbc.copy())

    @property
    def dtype(self):
        return self.array.dtype

    @property
    def size(self):
        return self.array.size

    @property
    def T(self):
        return self.array.T

    @property
    def flat(self):
        return self.array.flat

    @property
    def rank(self):
        """"Return the dimension of the cell.

        Equal to the number of nonzero lattice vectors."""
        # The name ndim clashes with ndarray.ndim
        return self.array.any(1).sum()

    @property
    def orthorhombic(self):
        """Return whether this cell is represented by a diagonal matrix."""
        from ase.geometry.cell import is_orthorhombic
        return is_orthorhombic(self.array)

    def lengths(self):
        """Return the length of each lattice vector as an array."""
        return np.array([np.linalg.norm(v) for v in self.array])

    def angles(self):
        """Return an array with the three angles alpha, beta, and gamma."""
        return self.cellpar()[3:].copy()

    @property
    def ndim(self):
        return self.array.ndim

    def __array__(self, dtype=float):
        if dtype != float:
            raise ValueError('Cannot convert cell to array of type {}'
                             .format(dtype))
        return self.array

    def __bool__(self):
        return bool(self.array.any())

    def __ne__(self, other):
        return self.array != other

    def __eq__(self, other):
        return self.array == other

    __nonzero__ = __bool__

    @property
    def volume(self):
        """Get the volume of this cell.

        If there are less than 3 lattice vectors, return 0."""
        # Fail or 0 for <3D cells?
        # Definitely 0 since this is currently a property.
        # I think normally it is more convenient just to get zero
        return np.abs(np.linalg.det(self.array))

    def tolist(self):
        return self.array.tolist()

    def scaled_positions(self, positions):
        """Calculate scaled positions from Cartesian positions."""
        return np.linalg.solve(self.complete().array.T, positions.T).T

    def cartesian_positions(self, scaled_positions):
        """Calculate Cartesian positions from scaled positions."""
        return scaled_positions @ self.complete()

    def reciprocal(self):
        """Get reciprocal lattice as a 3x3 array.

        Does not include factor of 2 pi."""
        return np.linalg.pinv(self.array).transpose()

    def __repr__(self):
        if self.orthorhombic:
            numbers = self.lengths().tolist()
        else:
            numbers = self.tolist()

        pbc = self.pbc
        if all(pbc):
            pbc = True
        elif not any(pbc):
            pbc = False
        return 'Cell({}, pbc={})'.format(numbers, pbc)

    def niggli_reduce(self, eps=1e-5):
        """Niggli reduce this cell, returning a new cell and mapping.

        See also :func:`ase.build.tools.niggli_reduce_cell`."""
        from ase.build.tools import niggli_reduce_cell
        cell, op = niggli_reduce_cell(self.array, epsfactor=1e-5)
        return Cell(cell, self.pbc), op

    def minkowski_reduce(self):
        """Minkowski-reduce this cell, returning new cell and mapping.

        See also :func:`ase.geometry.minkowski_reduction.minkowski_reduce`."""
        from ase.geometry.minkowski_reduction import minkowski_reduce
        rcell, op = minkowski_reduce(self)
        return Cell(rcell, self.pbc), op

    # XXX We want a reduction function that brings the cell into
    # standard form as defined by Setyawan and Curtarolo.
