# Copyright 2008, 2009 CAMd
# (see accompanying license files for details).

"""Definition of the Atoms class.

This module defines the central object in the ASE package: the Atoms
object.
"""

import copy
import numbers
import warnings
from math import cos, pi, sin

import numpy as np

import ase_koopmans.units as units
from ase_koopmans.atom import Atom
from ase_koopmans.cell import Cell
from ase_koopmans.constraints import (FixBondLengths, FixConstraint,
                                      FixLinearTriatomic)
from ase_koopmans.data import atomic_masses, atomic_masses_common
from ase_koopmans.geometry import (find_mic, get_angles, get_distances,
                                   wrap_positions)
from ase_koopmans.symbols import Symbols, symbols2numbers
from ase_koopmans.utils import deprecated


class Atoms(object):
    """Atoms object.

    The Atoms object can represent an isolated molecule, or a
    periodically repeated structure.  It has a unit cell and
    there may be periodic boundary conditions along any of the three
    unit cell axes.
    Information about the atoms (atomic numbers and position) is
    stored in ndarrays.  Optionally, there can be information about
    tags, momenta, masses, magnetic moments and charges.

    In order to calculate energies, forces and stresses, a calculator
    object has to attached to the atoms object.

    Parameters:

    symbols: str (formula) or list of str
        Can be a string formula, a list of symbols or a list of
        Atom objects.  Examples: 'H2O', 'COPt12', ['H', 'H', 'O'],
        [Atom('Ne', (x, y, z)), ...].
    positions: list of xyz-positions
        Atomic positions.  Anything that can be converted to an
        ndarray of shape (n, 3) will do: [(x1,y1,z1), (x2,y2,z2),
        ...].
    scaled_positions: list of scaled-positions
        Like positions, but given in units of the unit cell.
        Can not be set at the same time as positions.
    numbers: list of int
        Atomic numbers (use only one of symbols/numbers).
    tags: list of int
        Special purpose tags.
    momenta: list of xyz-momenta
        Momenta for all atoms.
    masses: list of float
        Atomic masses in atomic units.
    magmoms: list of float or list of xyz-values
        Magnetic moments.  Can be either a single value for each atom
        for collinear calculations or three numbers for each atom for
        non-collinear calculations.
    charges: list of float
        Initial atomic charges.
    cell: 3x3 matrix or length 3 or 6 vector
        Unit cell vectors.  Can also be given as just three
        numbers for orthorhombic cells, or 6 numbers, where
        first three are lengths of unit cell vectors, and the
        other three are angles between them (in degrees), in following order:
        [len(a), len(b), len(c), angle(b,c), angle(a,c), angle(a,b)].
        First vector will lie in x-direction, second in xy-plane,
        and the third one in z-positive subspace.
        Default value: [0, 0, 0].
    celldisp: Vector
        Unit cell displacement vector. To visualize a displaced cell
        around the center of mass of a Systems of atoms. Default value
        = (0,0,0)
    pbc: one or three bool
        Periodic boundary conditions flags.  Examples: True,
        False, 0, 1, (1, 1, 0), (True, False, False).  Default
        value: False.
    constraint: constraint object(s)
        Used for applying one or more constraints during structure
        optimization.
    calculator: calculator object
        Used to attach a calculator for calculating energies and atomic
        forces.
    info: dict of key-value pairs
        Dictionary of key-value pairs with additional information
        about the system.  The following keys may be used by ase_koopmans:

          - spacegroup: Spacegroup instance
          - unit_cell: 'conventional' | 'primitive' | int | 3 ints
          - adsorbate_info: Information about special adsorption sites

        Items in the info attribute survives copy and slicing and can
        be stored in and retrieved from trajectory files given that the
        key is a string, the value is JSON-compatible and, if the value is a
        user-defined object, its base_koopmans class is importable.  One should
        not make any assumptions about the existence of keys.

    Examples:

    These three are equivalent:

    >>> d = 1.104  # N2 bondlength
    >>> a = Atoms('N2', [(0, 0, 0), (0, 0, d)])
    >>> a = Atoms(numbers=[7, 7], positions=[(0, 0, 0), (0, 0, d)])
    >>> a = Atoms([Atom('N', (0, 0, 0)), Atom('N', (0, 0, d))])

    FCC gold:

    >>> a = 4.05  # Gold lattice constant
    >>> b = a / 2
    >>> fcc = Atoms('Au',
    ...             cell=[(0, b, b), (b, 0, b), (b, b, 0)],
    ...             pbc=True)

    Hydrogen wire:

    >>> d = 0.9  # H-H distance
    >>> h = Atoms('H', positions=[(0, 0, 0)],
    ...           cell=(d, 0, 0),
    ...           pbc=(1, 0, 0))
    """

    ase_koopmans_objtype = 'atoms'  # For JSONability

    def __init__(self, symbols=None,
                 positions=None, numbers=None,
                 tags=None, momenta=None, masses=None,
                 magmoms=None, charges=None,
                 scaled_positions=None,
                 cell=None, pbc=None, celldisp=None,
                 constraint=None,
                 calculator=None,
                 info=None,
                 velocities=None):

        self._cellobj = Cell.new()
        self._pbc = np.zeros(3, bool)

        atoms = None

        if hasattr(symbols, 'get_positions'):
            atoms = symbols
            symbols = None
        elif (isinstance(symbols, (list, tuple)) and
              len(symbols) > 0 and isinstance(symbols[0], Atom)):
            # Get data from a list or tuple of Atom objects:
            data = [[atom.get_raw(name) for atom in symbols]
                    for name in
                    ['position', 'number', 'tag', 'momentum',
                     'mass', 'magmom', 'charge']]
            atoms = self.__class__(None, *data)
            symbols = None

        if atoms is not None:
            # Get data from another Atoms object:
            if scaled_positions is not None:
                raise NotImplementedError
            if symbols is None and numbers is None:
                numbers = atoms.get_atomic_numbers()
            if positions is None:
                positions = atoms.get_positions()
            if tags is None and atoms.has('tags'):
                tags = atoms.get_tags()
            if momenta is None and atoms.has('momenta'):
                momenta = atoms.get_momenta()
            if magmoms is None and atoms.has('initial_magmoms'):
                magmoms = atoms.get_initial_magnetic_moments()
            if masses is None and atoms.has('masses'):
                masses = atoms.get_masses()
            if charges is None and atoms.has('initial_charges'):
                charges = atoms.get_initial_charges()
            if cell is None:
                cell = atoms.get_cell()
            if celldisp is None:
                celldisp = atoms.get_celldisp()
            if pbc is None:
                pbc = atoms.get_pbc()
            if constraint is None:
                constraint = [c.copy() for c in atoms.constraints]
            if calculator is None:
                calculator = atoms.calc
            if info is None:
                info = copy.deepcopy(atoms.info)

        self.arrays = {}

        if symbols is None:
            if numbers is None:
                if positions is not None:
                    natoms = len(positions)
                elif scaled_positions is not None:
                    natoms = len(scaled_positions)
                else:
                    natoms = 0
                numbers = np.zeros(natoms, int)
            self.new_array('numbers', numbers, int)
        else:
            if numbers is not None:
                raise TypeError(
                    'Use only one of "symbols" and "numbers".')
            else:
                self.new_array('numbers', symbols2numbers(symbols), int)

        if self.numbers.ndim != 1:
            raise ValueError('"numbers" must be 1-dimensional.')

        if cell is None:
            cell = np.zeros((3, 3))
        self.set_cell(cell)

        if celldisp is None:
            celldisp = np.zeros(shape=(3, 1))
        self.set_celldisp(celldisp)

        if positions is None:
            if scaled_positions is None:
                positions = np.zeros((len(self.arrays['numbers']), 3))
            else:
                assert self.number_of_lattice_vectors == 3
                positions = np.dot(scaled_positions, self.cell)
        else:
            if scaled_positions is not None:
                raise TypeError(
                    'Use only one of "symbols" and "numbers".')
        self.new_array('positions', positions, float, (3,))

        self.set_constraint(constraint)
        self.set_tags(default(tags, 0))
        self.set_masses(default(masses, None))
        self.set_initial_magnetic_moments(default(magmoms, 0.0))
        self.set_initial_charges(default(charges, 0.0))
        if pbc is None:
            pbc = False
        self.set_pbc(pbc)
        self.set_momenta(default(momenta, (0.0, 0.0, 0.0)),
                         apply_constraint=False)

        #                          V-- if instantiaed from list of Atom objs
        if velocities is not None and None not in velocities:
            if momenta is None:
                self.set_velocities(velocities)
            else:
                raise TypeError(
                    'Use only one of "momenta" and "velocities".')

        if info is None:
            self.info = {}
        else:
            self.info = dict(info)

        self.calc = calculator

    @property
    def symbols(self):
        """Get chemical symbols as a :class:`ase_koopmans.symbols.Symbols` object.

        The object works like ``atoms.numbers`` except its values
        are strings.  It supports in-place editing."""
        return Symbols(self.numbers)

    @symbols.setter
    def symbols(self, obj):
        new_symbols = Symbols.fromsymbols(obj)
        self.numbers[:] = new_symbols.numbers

    @deprecated(DeprecationWarning('Please_koopmans use atoms.calc = calc'))
    def set_calculator(self, calc=None):
        """Attach calculator object.

        Please_koopmans use the equivalent atoms.calc = calc instead of this
        method."""
        self.calc = calc

    @deprecated(DeprecationWarning('Please_koopmans use atoms.calc'))
    def get_calculator(self):
        """Get currently attached calculator object.

        Please_koopmans use the equivalent atoms.calc instead of
        atoms.get_calculator()."""
        return self.calc

    @property
    def calc(self):
        """Calculator object."""
        return self._calc

    @calc.setter
    def calc(self, calc):
        self._calc = calc
        if hasattr(calc, 'set_atoms'):
            calc.set_atoms(self)

    @calc.deleter  # type: ignore
    @deprecated(DeprecationWarning('Please_koopmans use atoms.calc = None'))
    def calc(self):
        self._calc = None

    @property
    def number_of_lattice_vectors(self):
        """Number of (non-zero) lattice vectors."""
        return self.cell.rank

    def set_constraint(self, constraint=None):
        """Apply one or more constrains.

        The *constraint* argument must be one constraint object or a
        list of constraint objects."""
        if constraint is None:
            self._constraints = []
        else:
            if isinstance(constraint, list):
                self._constraints = constraint
            elif isinstance(constraint, tuple):
                self._constraints = list(constraint)
            else:
                self._constraints = [constraint]

    def _get_constraints(self):
        return self._constraints

    def _del_constraints(self):
        self._constraints = []

    constraints = property(_get_constraints, set_constraint, _del_constraints,
                           'Constraints of the atoms.')

    def set_cell(self, cell, scale_atoms=False, apply_constraint=True):
        """Set unit cell vectors.

        Parameters:

        cell: 3x3 matrix or length 3 or 6 vector
            Unit cell.  A 3x3 matrix (the three unit cell vectors) or
            just three numbers for an orthorhombic cell. Another option is
            6 numbers, which describes unit cell with lengths of unit cell
            vectors and with angles between them (in degrees), in following
            order: [len(a), len(b), len(c), angle(b,c), angle(a,c),
            angle(a,b)].  First vector will lie in x-direction, second in
            xy-plane, and the third one in z-positive subspace.
        scale_atoms: bool
            Fix atomic positions or move atoms with the unit cell?
            Default behavior is to *not* move the atoms (scale_atoms=False).

        Examples:

        Two equivalent ways to define an orthorhombic cell:

        >>> atoms = Atoms('He')
        >>> a, b, c = 7, 7.5, 8
        >>> atoms.set_cell([a, b, c])
        >>> atoms.set_cell([(a, 0, 0), (0, b, 0), (0, 0, c)])

        FCC unit cell:

        >>> atoms.set_cell([(0, b, b), (b, 0, b), (b, b, 0)])

        Hexagonal unit cell:

        >>> atoms.set_cell([a, a, c, 90, 90, 120])

        Rhombohedral unit cell:

        >>> alpha = 77
        >>> atoms.set_cell([a, a, a, alpha, alpha, alpha])
        """

        # Override pbcs if and only if given a Cell object:
        cell = Cell.new(cell)

        # XXX not working well during initialize due to missing _constraints
        if apply_constraint and hasattr(self, '_constraints'):
            for constraint in self.constraints:
                if hasattr(constraint, 'adjust_cell'):
                    constraint.adjust_cell(self, cell)

        if scale_atoms:
            M = np.linalg.solve(self.cell.complete(), cell.complete())
            self.positions[:] = np.dot(self.positions, M)

        self.cell[:] = cell

    def set_celldisp(self, celldisp):
        """Set the unit cell displacement vectors."""
        celldisp = np.array(celldisp, float)
        self._celldisp = celldisp

    def get_celldisp(self):
        """Get the unit cell displacement vectors."""
        return self._celldisp.copy()

    def get_cell(self, complete=False):
        """Get the three unit cell vectors as a `class`:ase_koopmans.cell.Cell` object.

        The Cell object resembles a 3x3 ndarray, and cell[i, j]
        is the jth Cartesian coordinate of the ith cell vector."""
        if complete:
            cell = self.cell.complete()
        else:
            cell = self.cell.copy()

        return cell

    def get_cell_lengths_and_angles(self):
        """Get unit cell parameters. Sequence of 6 numbers.

        First three are unit cell vector lengths and second three
        are angles between them::

            [len(a), len(b), len(c), angle(b,c), angle(a,c), angle(a,b)]

        in degrees.
        """
        return self.cell.cellpar()

    def get_reciprocal_cell(self):
        """Get the three reciprocal lattice vectors as a 3x3 ndarray.

        Note that the commonly used factor of 2 pi for Fourier
        transforms is not included here."""

        return self.cell.reciprocal()

    def set_pbc(self, pbc):
        """Set periodic boundary condition flags."""
        self._pbc[:] = pbc

    def get_pbc(self):
        """Get periodic boundary condition flags."""
        return self.pbc.copy()

    def new_array(self, name, a, dtype=None, shape=None):
        """Add new array.

        If *shape* is not *None*, the shape of *a* will be checked."""

        if dtype is not None:
            a = np.array(a, dtype, order='C')
            if len(a) == 0 and shape is not None:
                a.shape = (-1,) + shape
        else:
            if not a.flags['C_CONTIGUOUS']:
                a = np.ascontiguousarray(a)
            else:
                a = a.copy()

        if name in self.arrays:
            raise RuntimeError('Array {} already present'.format(name))

        for b in self.arrays.values():
            if len(a) != len(b):
                raise ValueError('Array "%s" has wrong length: %d != %d.' %
                                 (name, len(a), len(b)))
            break

        if shape is not None and a.shape[1:] != shape:
            raise ValueError('Array "%s" has wrong shape %s != %s.' %
                             (a.name, a.shape, (a.shape[0:1] + shape)))

        self.arrays[name] = a

    def get_array(self, name, copy=True):
        """Get an array.

        Returns a copy unless the optional argument copy is false.
        """
        if copy:
            return self.arrays[name].copy()
        else:
            return self.arrays[name]

    def set_array(self, name, a, dtype=None, shape=None):
        """Update array.

        If *shape* is not *None*, the shape of *a* will be checked.
        If *a* is *None*, then the array is deleted."""

        b = self.arrays.get(name)
        if b is None:
            if a is not None:
                self.new_array(name, a, dtype, shape)
        else:
            if a is None:
                del self.arrays[name]
            else:
                a = np.asarray(a)
                if a.shape != b.shape:
                    raise ValueError('Array "%s" has wrong shape %s != %s.' %
                                     (name, a.shape, b.shape))
                b[:] = a

    def has(self, name):
        """Check for existence of array.

        name must be one of: 'tags', 'momenta', 'masses', 'initial_magmoms',
        'initial_charges'."""
        # XXX extend has to calculator properties
        return name in self.arrays

    def set_atomic_numbers(self, numbers):
        """Set atomic numbers."""
        self.set_array('numbers', numbers, int, ())

    def get_atomic_numbers(self):
        """Get integer array of atomic numbers."""
        return self.arrays['numbers'].copy()

    def get_chemical_symbols(self):
        """Get list of chemical symbol strings.

        Equivalent to ``list(atoms.symbols)``."""
        return list(self.symbols)

    def set_chemical_symbols(self, symbols):
        """Set chemical symbols."""
        self.set_array('numbers', symbols2numbers(symbols), int, ())

    def get_chemical_formula(self, mode='hill', empirical=False):
        """Get the chemical formula as a string base_koopmansd on the chemical symbols.

        Parameters:

        mode: str
            There are four different modes available:

            'all': The list of chemical symbols are contracted to a string,
            e.g. ['C', 'H', 'H', 'H', 'O', 'H'] becomes 'CHHHOH'.

            'reduce': The same as 'all' where repeated elements are contracted
            to a single symbol and a number, e.g. 'CHHHOCHHH' is reduced to
            'CH3OCH3'.

            'hill': The list of chemical symbols are contracted to a string
            following the Hill notation (alphabetical order with C and H
            first), e.g. 'CHHHOCHHH' is reduced to 'C2H6O' and 'SOOHOHO' to
            'H2O4S'. This is default.

            'metal': The list of chemical symbols (alphabetical metals,
            and alphabetical non-metals)

        empirical, bool (optional, default=False)
            Divide the symbol counts by their greatest common divisor to yield
            an empirical formula. Only for mode `metal` and `hill`.
        """
        return self.symbols.get_chemical_formula(mode, empirical)

    def set_tags(self, tags):
        """Set tags for all atoms. If only one tag is supplied, it is
        applied to all atoms."""
        if isinstance(tags, int):
            tags = [tags] * len(self)
        self.set_array('tags', tags, int, ())

    def get_tags(self):
        """Get integer array of tags."""
        if 'tags' in self.arrays:
            return self.arrays['tags'].copy()
        else:
            return np.zeros(len(self), int)

    def set_momenta(self, momenta, apply_constraint=True):
        """Set momenta."""
        if (apply_constraint and len(self.constraints) > 0 and
           momenta is not None):
            momenta = np.array(momenta)  # modify a copy
            for constraint in self.constraints:
                if hasattr(constraint, 'adjust_momenta'):
                    constraint.adjust_momenta(self, momenta)
        self.set_array('momenta', momenta, float, (3,))

    def set_velocities(self, velocities):
        """Set the momenta by specifying the velocities."""
        self.set_momenta(self.get_masses()[:, np.newaxis] * velocities)

    def get_momenta(self):
        """Get array of momenta."""
        if 'momenta' in self.arrays:
            return self.arrays['momenta'].copy()
        else:
            return np.zeros((len(self), 3))

    def set_masses(self, masses='defaults'):
        """Set atomic masses in atomic mass units.

        The array masses should contain a list of masses.  In case_koopmans
        the masses argument is not given or for those elements of the
        masses list that are None, standard values are set."""

        if isinstance(masses, str):
            if masses == 'defaults':
                masses = atomic_masses[self.arrays['numbers']]
            elif masses == 'most_common':
                masses = atomic_masses_common[self.arrays['numbers']]
        elif isinstance(masses, (list, tuple)):
            newmasses = []
            for m, Z in zip(masses, self.arrays['numbers']):
                if m is None:
                    newmasses.append(atomic_masses[Z])
                else:
                    newmasses.append(m)
            masses = newmasses
        self.set_array('masses', masses, float, ())

    def get_masses(self):
        """Get array of masses in atomic mass units."""
        if 'masses' in self.arrays:
            return self.arrays['masses'].copy()
        else:
            return atomic_masses[self.arrays['numbers']]

    def set_initial_magnetic_moments(self, magmoms=None):
        """Set the initial magnetic moments.

        Use either one or three numbers for every atom (collinear
        or non-collinear spins)."""

        if magmoms is None:
            self.set_array('initial_magmoms', None)
        else:
            magmoms = np.asarray(magmoms)
            self.set_array('initial_magmoms', magmoms, float,
                           magmoms.shape[1:])

    def get_initial_magnetic_moments(self):
        """Get array of initial magnetic moments."""
        if 'initial_magmoms' in self.arrays:
            return self.arrays['initial_magmoms'].copy()
        else:
            return np.zeros(len(self))

    def get_magnetic_moments(self):
        """Get calculated local magnetic moments."""
        if self._calc is None:
            raise RuntimeError('Atoms object has no calculator.')
        return self._calc.get_magnetic_moments(self)

    def get_magnetic_moment(self):
        """Get calculated total magnetic moment."""
        if self._calc is None:
            raise RuntimeError('Atoms object has no calculator.')
        return self._calc.get_magnetic_moment(self)

    def set_initial_charges(self, charges=None):
        """Set the initial charges."""

        if charges is None:
            self.set_array('initial_charges', None)
        else:
            self.set_array('initial_charges', charges, float, ())

    def get_initial_charges(self):
        """Get array of initial charges."""
        if 'initial_charges' in self.arrays:
            return self.arrays['initial_charges'].copy()
        else:
            return np.zeros(len(self))

    def get_charges(self):
        """Get calculated charges."""
        if self._calc is None:
            raise RuntimeError('Atoms object has no calculator.')
        try:
            return self._calc.get_charges(self)
        except AttributeError:
            from ase_koopmans.calculators.calculator import \
                PropertyNotImplementedError
            raise PropertyNotImplementedError

    def set_positions(self, newpositions, apply_constraint=True):
        """Set positions, honoring any constraints. To ignore constraints,
        use *apply_constraint=False*."""
        if self.constraints and apply_constraint:
            newpositions = np.array(newpositions, float)
            for constraint in self.constraints:
                constraint.adjust_positions(self, newpositions)

        self.set_array('positions', newpositions, shape=(3,))

    def get_positions(self, wrap=False, **wrap_kw):
        """Get array of positions.

        Parameters:

        wrap: bool
            wrap atoms back to the cell before returning positions
        wrap_kw: (keyword=value) pairs
            optional keywords `pbc`, `center`, `pretty_translation`, `eps`,
            see :func:`ase_koopmans.geometry.wrap_positions`
        """
        if wrap:
            if 'pbc' not in wrap_kw:
                wrap_kw['pbc'] = self.pbc
            return wrap_positions(self.positions, self.cell, **wrap_kw)
        else:
            return self.arrays['positions'].copy()

    def get_potential_energy(self, force_consistent=False,
                             apply_constraint=True):
        """Calculate potential energy.

        Ask the attached calculator to calculate the potential energy and
        apply constraints.  Use *apply_constraint=False* to get the raw
        forces.

        When supported by the calculator, either the energy extrapolated
        to zero Kelvin or the energy consistent with the forces (the free
        energy) can be returned.
        """
        if self._calc is None:
            raise RuntimeError('Atoms object has no calculator.')
        if force_consistent:
            energy = self._calc.get_potential_energy(
                self, force_consistent=force_consistent)
        else:
            energy = self._calc.get_potential_energy(self)
        if apply_constraint:
            for constraint in self.constraints:
                if hasattr(constraint, 'adjust_potential_energy'):
                    energy += constraint.adjust_potential_energy(self)
        return energy

    def get_properties(self, properties):
        """This method is experimental; currently for internal use."""
        # XXX Something about constraints.
        if self._calc is None:
            raise RuntimeError('Atoms object has no calculator.')
        return self._calc.calculate_properties(self, properties)

    def get_potential_energies(self):
        """Calculate the potential energies of all the atoms.

        Only available with calculators supporting per-atom energies
        (e.g. classical potentials).
        """
        if self._calc is None:
            raise RuntimeError('Atoms object has no calculator.')
        return self._calc.get_potential_energies(self)

    def get_kinetic_energy(self):
        """Get the kinetic energy."""
        momenta = self.arrays.get('momenta')
        if momenta is None:
            return 0.0
        return 0.5 * np.vdot(momenta, self.get_velocities())

    def get_velocities(self):
        """Get array of velocities."""
        momenta = self.arrays.get('momenta')
        if momenta is None:
            return None
        m = self.arrays.get('masses')
        if m is None:
            m = atomic_masses[self.arrays['numbers']]
        return momenta / m.reshape(-1, 1)

    def get_total_energy(self):
        """Get the total energy - potential plus kinetic energy."""
        return self.get_potential_energy() + self.get_kinetic_energy()

    def get_forces(self, apply_constraint=True, md=False):
        """Calculate atomic forces.

        Ask the attached calculator to calculate the forces and apply
        constraints.  Use *apply_constraint=False* to get the raw
        forces.

        For molecular dynamics (md=True) we don't apply the constraint
        to the forces but to the momenta. When holonomic constraints for
        rigid linear triatomic molecules are present, ask the constraints
        to redistribute the forces within each triple defined in the
        constraints (required for molecular dynamics with this type of
        constraints)."""

        if self._calc is None:
            raise RuntimeError('Atoms object has no calculator.')
        forces = self._calc.get_forces(self)

        if apply_constraint:
            # We need a special md flag here because for MD we want
            # to skip real constraints but include special "constraints"
            # Like Hookean.
            for constraint in self.constraints:
                if md and hasattr(constraint, 'redistribute_forces_md'):
                    constraint.redistribute_forces_md(self, forces)
                if not md or hasattr(constraint, 'adjust_potential_energy'):
                    constraint.adjust_forces(self, forces)
        return forces

    # Informs calculators (e.g. Asap) that ideal gas contribution is added here.
    _ase_koopmans_handles_dynamic_stress = True

    def get_stress(self, voigt=True, apply_constraint=True,
                   include_ideal_gas=False):
        """Calculate stress tensor.

        Returns an array of the six independent components of the
        symmetric stress tensor, in the traditional Voigt order
        (xx, yy, zz, yz, xz, xy) or as a 3x3 matrix.  Default is Voigt
        order.

        The ideal gas contribution to the stresses is added if the
        atoms have momenta and ``include_ideal_gas`` is set to True.
        """

        if self._calc is None:
            raise RuntimeError('Atoms object has no calculator.')

        stress = self._calc.get_stress(self)
        shape = stress.shape

        if shape == (3, 3):
            # Convert to the Voigt form before possibly applying
            # constraints and adding the dynamic part of the stress
            # (the "ideal gas contribution").
            stress = np.array([stress[0, 0], stress[1, 1], stress[2, 2],
                               stress[1, 2], stress[0, 2], stress[0, 1]])
        else:
            assert shape == (6,)

        if apply_constraint:
            for constraint in self.constraints:
                if hasattr(constraint, 'adjust_stress'):
                    constraint.adjust_stress(self, stress)

        # Add ideal gas contribution, if applicable
        if include_ideal_gas and self.has('momenta'):
            stresscomp = np.array([[0, 5, 4], [5, 1, 3], [4, 3, 2]])
            p = self.get_momenta()
            masses = self.get_masses()
            invmass = 1.0 / masses
            invvol = 1.0 / self.get_volume()
            for alpha in range(3):
                for beta in range(alpha, 3):
                    stress[stresscomp[alpha, beta]] -= (
                        p[:, alpha] * p[:, beta] * invmass).sum() * invvol

        if voigt:
            return stress
        else:
            xx, yy, zz, yz, xz, xy = stress
            return np.array([(xx, xy, xz),
                             (xy, yy, yz),
                             (xz, yz, zz)])

    def get_stresses(self, include_ideal_gas=False):
        """Calculate the stress-tensor of all the atoms.

        Only available with calculators supporting per-atom energies and
        stresses (e.g. classical potentials).  Even for such calculators
        there is a certain arbitrariness in defining per-atom stresses.

        The ideal gas contribution to the stresses is added if the
        atoms have momenta and ``include_ideal_gas`` is set to True.
        """
        if self._calc is None:
            raise RuntimeError('Atoms object has no calculator.')
        stresses = self._calc.get_stresses(self)
        if include_ideal_gas and self.has('momenta'):
            stresscomp = np.array([[0, 5, 4], [5, 1, 3], [4, 3, 2]])
            if hasattr(self._calc, 'get_atomic_volumes'):
                invvol = 1.0 / self._calc.get_atomic_volumes()
            else:
                invvol = self.get_global_number_of_atoms() / self.get_volume()
            p = self.get_momenta()
            invmass = 1.0 / self.get_masses()
            for alpha in range(3):
                for beta in range(alpha, 3):
                    stresses[:, stresscomp[alpha, beta]] -= (
                        p[:, alpha] * p[:, beta] * invmass * invvol)
        return stresses

    def get_dipole_moment(self):
        """Calculate the electric dipole moment for the atoms object.

        Only available for calculators which has a get_dipole_moment()
        method."""

        if self._calc is None:
            raise RuntimeError('Atoms object has no calculator.')
        return self._calc.get_dipole_moment(self)

    def copy(self):
        """Return a copy."""
        atoms = self.__class__(cell=self.cell, pbc=self.pbc, info=self.info,
                               celldisp=self._celldisp.copy())

        atoms.arrays = {}
        for name, a in self.arrays.items():
            atoms.arrays[name] = a.copy()
        atoms.constraints = copy.deepcopy(self.constraints)
        return atoms

    def todict(self):
        """For basic JSON (non-database_koopmans) support."""
        d = dict(self.arrays)
        d['cell'] = np.asarray(self.cell)
        d['pbc'] = self.pbc
        if self._celldisp.any():
            d['celldisp'] = self._celldisp
        if self.constraints:
            d['constraints'] = self.constraints
        if self.info:
            d['info'] = self.info
        # Calculator...  trouble.
        return d

    @classmethod
    def fromdict(cls, dct):
        """Rebuild atoms object from dictionary representation (todict)."""
        dct = dct.copy()
        kw = {}
        for name in ['numbers', 'positions', 'cell', 'pbc']:
            kw[name] = dct.pop(name)

        constraints = dct.pop('constraints', None)
        if constraints:
            from ase_koopmans.constraints import dict2constraint
            constraints = [dict2constraint(d) for d in constraints]

        atoms = cls(constraint=constraints,
                    celldisp=dct.pop('celldisp', None),
                    info=dct.pop('info', None), **kw)
        natoms = len(atoms)

        # Some arrays are named differently from the atoms __init__ keywords.
        # Also, there may be custom arrays.  Hence we set them directly:
        for name, arr in dct.items():
            assert len(arr) == natoms, name
            assert isinstance(arr, np.ndarray)
            atoms.arrays[name] = arr
        return atoms

    def __len__(self):
        return len(self.arrays['positions'])

    def get_number_of_atoms(self):
        """Deprecated, please_koopmans do not use.

        You probably want len(atoms).  Or if your atoms are distributed,
        use (and see) get_global_number_of_atoms()."""
        import warnings
        warnings.warn('Use get_global_number_of_atoms() instead',
                      np.VisibleDeprecationWarning)
        return len(self)

    def get_global_number_of_atoms(self):
        """Returns the global number of atoms in a distributed-atoms parallel
        simulation.

        DO NOT USE UNLESS YOU KNOW WHAT YOU ARE DOING!

        Equivalent to len(atoms) in the standard ASE Atoms class.  You should
        normally use len(atoms) instead.  This function's only purpose is to
        make compatibility between ASE and Asap easier to maintain by having a
        few places in ASE use this function instead.  It is typically only
        when counting the global number of degrees of freedom or in similar
        situations.
        """
        return len(self)

    def __repr__(self):
        tokens = []

        N = len(self)
        if N <= 60:
            symbols = self.get_chemical_formula('reduce')
        else:
            symbols = self.get_chemical_formula('hill')
        tokens.append("symbols='{0}'".format(symbols))

        if self.pbc.any() and not self.pbc.all():
            tokens.append('pbc={0}'.format(self.pbc.tolist()))
        else:
            tokens.append('pbc={0}'.format(self.pbc[0]))

        cell = self.cell
        if cell:
            if cell.orthorhombic:
                cell = cell.lengths().tolist()
            else:
                cell = cell.tolist()
            tokens.append('cell={0}'.format(cell))

        for name in sorted(self.arrays):
            if name in ['numbers', 'positions']:
                continue
            tokens.append('{0}=...'.format(name))

        if self.constraints:
            if len(self.constraints) == 1:
                constraint = self.constraints[0]
            else:
                constraint = self.constraints
            tokens.append('constraint={0}'.format(repr(constraint)))

        if self._calc is not None:
            tokens.append('calculator={0}(...)'
                          .format(self._calc.__class__.__name__))

        return '{0}({1})'.format(self.__class__.__name__, ', '.join(tokens))

    def __add__(self, other):
        atoms = self.copy()
        atoms += other
        return atoms

    def extend(self, other):
        """Extend atoms object by appending atoms from *other*."""
        if isinstance(other, Atom):
            other = self.__class__([other])

        n1 = len(self)
        n2 = len(other)

        for name, a1 in self.arrays.items():
            a = np.zeros((n1 + n2,) + a1.shape[1:], a1.dtype)
            a[:n1] = a1
            if name == 'masses':
                a2 = other.get_masses()
            else:
                a2 = other.arrays.get(name)
            if a2 is not None:
                a[n1:] = a2
            self.arrays[name] = a

        for name, a2 in other.arrays.items():
            if name in self.arrays:
                continue
            a = np.empty((n1 + n2,) + a2.shape[1:], a2.dtype)
            a[n1:] = a2
            if name == 'masses':
                a[:n1] = self.get_masses()[:n1]
            else:
                a[:n1] = 0

            self.set_array(name, a)

    def __iadd__(self, other):
        self.extend(other)
        return self

    def append(self, atom):
        """Append atom to end."""
        self.extend(self.__class__([atom]))

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]

    def __getitem__(self, i):
        """Return a subset of the atoms.

        i -- scalar integer, list of integers, or slice object
        describing which atoms to return.

        If i is a scalar, return an Atom object. If i is a list or a
        slice, return an Atoms object with the same cell, pbc, and
        other associated info as the original Atoms object. The
        indices of the constraints will be shuffled so that they match
        the indexing in the subset returned.

        """

        if isinstance(i, numbers.Integral):
            natoms = len(self)
            if i < -natoms or i >= natoms:
                raise IndexError('Index out of range.')

            return Atom(atoms=self, index=i)
        elif not isinstance(i, slice):
            i = np.array(i)
            # if i is a mask
            if i.dtype == bool:
                if len(i) != len(self):
                    raise IndexError('Length of mask {} must equal '
                                     'number of atoms {}'
                                     .format(len(i), len(self)))
                i = np.arange(len(self))[i]

        import copy

        conadd = []
        # Constraints need to be deepcopied, but only the relevant ones.
        for con in copy.deepcopy(self.constraints):
            if isinstance(con, (FixConstraint, FixBondLengths,
                                FixLinearTriatomic)):
                try:
                    con.index_shuffle(self, i)
                    conadd.append(con)
                except IndexError:
                    pass

        atoms = self.__class__(cell=self.cell, pbc=self.pbc, info=self.info,
                               # should be communicated to the slice as well
                               celldisp=self._celldisp)
        # TODO: Do we need to shuffle indices in adsorbate_info too?

        atoms.arrays = {}
        for name, a in self.arrays.items():
            atoms.arrays[name] = a[i].copy()

        atoms.constraints = conadd
        return atoms

    def __delitem__(self, i):
        from ase_koopmans.constraints import FixAtoms
        for c in self._constraints:
            if not isinstance(c, FixAtoms):
                raise RuntimeError('Remove constraint using set_constraint() '
                                   'before deleting atoms.')

        if isinstance(i, list) and len(i) > 0:
            # Make sure a list of booleans will work correctly and not be
            # interpreted at 0 and 1 indices.
            i = np.array(i)

        if len(self._constraints) > 0:
            n = len(self)
            i = np.arange(n)[i]
            if isinstance(i, int):
                i = [i]
            constraints = []
            for c in self._constraints:
                c = c.delete_atoms(i, n)
                if c is not None:
                    constraints.append(c)
            self.constraints = constraints

        mask = np.ones(len(self), bool)
        mask[i] = False
        for name, a in self.arrays.items():
            self.arrays[name] = a[mask]

    def pop(self, i=-1):
        """Remove and return atom at index *i* (default last)."""
        atom = self[i]
        atom.cut_reference_to_atoms()
        del self[i]
        return atom

    def __imul__(self, m):
        """In-place repeat of atoms."""
        if isinstance(m, int):
            m = (m, m, m)

        for x, vec in zip(m, self.cell):
            if x != 1 and not vec.any():
                raise ValueError('Cannot repeat along undefined lattice '
                                 'vector')

        M = np.prod(m)
        n = len(self)

        for name, a in self.arrays.items():
            self.arrays[name] = np.tile(a, (M,) + (1,) * (len(a.shape) - 1))

        positions = self.arrays['positions']
        i0 = 0
        for m0 in range(m[0]):
            for m1 in range(m[1]):
                for m2 in range(m[2]):
                    i1 = i0 + n
                    positions[i0:i1] += np.dot((m0, m1, m2), self.cell)
                    i0 = i1

        if self.constraints is not None:
            self.constraints = [c.repeat(m, n) for c in self.constraints]

        self.cell = np.array([m[c] * self.cell[c] for c in range(3)])

        return self

    def repeat(self, rep):
        """Create new repeated atoms object.

        The *rep* argument should be a sequence of three positive
        integers like *(2,3,1)* or a single integer (*r*) equivalent
        to *(r,r,r)*."""

        atoms = self.copy()
        atoms *= rep
        return atoms

    def __mul__(self, rep):
        return self.repeat(rep)

    def translate(self, displacement):
        """Translate atomic positions.

        The displacement argument can be a float an xyz vector or an
        nx3 array (where n is the number of atoms)."""

        self.arrays['positions'] += np.array(displacement)

    def center(self, vacuum=None, axis=(0, 1, 2), about=None):
        """Center atoms in unit cell.

        Centers the atoms in the unit cell, so there is the same
        amount of vacuum on all sides.

        vacuum: float (default: None)
            If specified adjust the amount of vacuum when centering.
            If vacuum=10.0 there will thus be 10 Angstrom of vacuum
            on each side.
        axis: int or sequence of ints
            Axis or axes to act on.  Default: Act on all axes.
        about: float or array (default: None)
            If specified, center the atoms about <about>.
            I.e., about=(0., 0., 0.) (or just "about=0.", interpreted
            identically), to center about the origin.
        """

        # Find the orientations of the faces of the unit cell
        cell = self.cell.complete()
        dirs = np.zeros_like(cell)
        for i in range(3):
            dirs[i] = np.cross(cell[i - 1], cell[i - 2])
            dirs[i] /= np.sqrt(np.dot(dirs[i], dirs[i]))  # normalize
            if np.dot(dirs[i], cell[i]) < 0.0:
                dirs[i] *= -1

        if isinstance(axis, int):
            axes = (axis,)
        else:
            axes = axis

        # if vacuum and any(self.pbc[x] for x in axes):
        #     warnings.warn(
        #         'You are adding vacuum along a periodic direction!')

        # Now, decide how much each basis vector should be made longer
        p = self.arrays['positions']
        longer = np.zeros(3)
        shift = np.zeros(3)
        for i in axes:
            p0 = np.dot(p, dirs[i]).min() if len(p) else 0
            p1 = np.dot(p, dirs[i]).max() if len(p) else 0
            height = np.dot(cell[i], dirs[i])
            if vacuum is not None:
                lng = (p1 - p0 + 2 * vacuum) - height
            else:
                lng = 0.0  # Do not change unit cell size!
            top = lng + height - p1
            shf = 0.5 * (top - p0)
            cosphi = np.dot(cell[i], dirs[i]) / np.sqrt(np.dot(cell[i],
                                                               cell[i]))
            longer[i] = lng / cosphi
            shift[i] = shf / cosphi

        # Now, do it!
        translation = np.zeros(3)
        for i in axes:
            nowlen = np.sqrt(np.dot(cell[i], cell[i]))
            if vacuum is not None or self.cell[i].any():
                self.cell[i] = cell[i] * (1 + longer[i] / nowlen)
                translation += shift[i] * cell[i] / nowlen
        self.arrays['positions'] += translation

        # Optionally, translate to center about a point in space.
        if about is not None:
            for vector in self.cell:
                self.positions -= vector / 2.0
            self.positions += about

    def get_center_of_mass(self, scaled=False):
        """Get the center of mass.

        If scaled=True the center of mass in scaled coordinates
        is returned."""
        m = self.get_masses()
        com = np.dot(m, self.arrays['positions']) / m.sum()
        if scaled:
            return np.linalg.solve(self.cell.T, com)
        else:
            return com

    def get_moments_of_inertia(self, vectors=False):
        """Get the moments of inertia along the principal axes.

        The three principal moments of inertia are computed from the
        eigenvalues of the symmetric inertial tensor. Periodic boundary
        conditions are ignored. Units of the moments of inertia are
        amu*angstrom**2.
        """
        com = self.get_center_of_mass()
        positions = self.get_positions()
        positions -= com  # translate center of mass to origin
        masses = self.get_masses()

        # Initialize elements of the inertial tensor
        I11 = I22 = I33 = I12 = I13 = I23 = 0.0
        for i in range(len(self)):
            x, y, z = positions[i]
            m = masses[i]

            I11 += m * (y ** 2 + z ** 2)
            I22 += m * (x ** 2 + z ** 2)
            I33 += m * (x ** 2 + y ** 2)
            I12 += -m * x * y
            I13 += -m * x * z
            I23 += -m * y * z

        I = np.array([[I11, I12, I13],
                      [I12, I22, I23],
                      [I13, I23, I33]])

        evals, evecs = np.linalg.eigh(I)
        if vectors:
            return evals, evecs.transpose()
        else:
            return evals

    def get_angular_momentum(self):
        """Get total angular momentum with respect to the center of mass."""
        com = self.get_center_of_mass()
        positions = self.get_positions()
        positions -= com  # translate center of mass to origin
        return np.cross(positions, self.get_momenta()).sum(0)

    def rotate(self, a, v=None, center=(0, 0, 0), rotate_cell=False):
        """Rotate atoms base_koopmansd on a vector and an angle, or two vectors.

        Parameters:

        a = None:
            Angle that the atoms is rotated around the vecor 'v'. 'a'
            can also be a vector and then 'a' is rotated
            into 'v'.

        v:
            Vector to rotate the atoms around. Vectors can be given as
            strings: 'x', '-x', 'y', ... .

        center = (0, 0, 0):
            The center is kept fixed under the rotation. Use 'COM' to fix
            the center of mass, 'COP' to fix the center of positions or
            'COU' to fix the center of cell.

        rotate_cell = False:
            If true the cell is also rotated.

        Examples:

        Rotate 90 degrees around the z-axis, so that the x-axis is
        rotated into the y-axis:

        >>> atoms = Atoms()
        >>> atoms.rotate(90, 'z')
        >>> atoms.rotate(90, (0, 0, 1))
        >>> atoms.rotate(-90, '-z')
        >>> atoms.rotate('x', 'y')
        >>> atoms.rotate((1, 0, 0), (0, 1, 0))
        """

        if not isinstance(a, (float, int)):
            # old API maybe?
            warning = ('Please_koopmans use new API: '
                       'atoms_obj.rotate(a, v) '
                       'where v is a vector to rotate around and '
                       'a is the angle in degrees.')
            if isinstance(v, (float, int)):
                warnings.warn(warning, FutureWarning)
                a, v = v * 180 / pi, a
            elif v is None:
                warnings.warn(warning, FutureWarning)
                v = a
                a = None
            else:
                assert a is not None
                a, v = v, a
        else:
            assert a is not None

        norm = np.linalg.norm
        v = string2vector(v)
        if a is None:
            a = norm(v) * 180 / pi  # old API

        normv = norm(v)

        if normv == 0.0:
            raise ZeroDivisionError('Cannot rotate: norm(v) == 0')

        if isinstance(a, (float, int)):
            a *= pi / 180
            v /= normv
            c = cos(a)
            s = sin(a)
        else:
            v2 = string2vector(a)
            v /= normv
            normv2 = np.linalg.norm(v2)
            if normv2 == 0:
                raise ZeroDivisionError('Cannot rotate: norm(a) == 0')
            v2 /= norm(v2)
            c = np.dot(v, v2)
            v = np.cross(v, v2)
            s = norm(v)
            # In case_koopmans *v* and *a* are parallel, np.cross(v, v2) vanish
            # and can't be used as a rotation axis. However, in this
            # case_koopmans any rotation axis perpendicular to v2 will do.
            eps = 1e-7
            if s < eps:
                v = np.cross((0, 0, 1), v2)
                if norm(v) < eps:
                    v = np.cross((1, 0, 0), v2)
                assert norm(v) >= eps
            elif s > 0:
                v /= s

        if isinstance(center, str):
            if center.lower() == 'com':
                center = self.get_center_of_mass()
            elif center.lower() == 'cop':
                center = self.get_positions().mean(axis=0)
            elif center.lower() == 'cou':
                center = self.get_cell().sum(axis=0) / 2
            else:
                raise ValueError('Cannot interpret center')
        else:
            center = np.array(center)

        p = self.arrays['positions'] - center
        self.arrays['positions'][:] = (c * p -
                                       np.cross(p, s * v) +
                                       np.outer(np.dot(p, v), (1.0 - c) * v) +
                                       center)
        if rotate_cell:
            rotcell = self.get_cell()
            rotcell[:] = (c * rotcell -
                          np.cross(rotcell, s * v) +
                          np.outer(np.dot(rotcell, v), (1.0 - c) * v))
            self.set_cell(rotcell)

    def rotate_euler(self, center=(0, 0, 0), phi=0.0, theta=0.0, psi=0.0):
        warnings.warn(
            'Please_koopmans use this method instead: '
            'euler_rotate(phi=0, theta=0, psi=0, center=(0, 0, 0)) '
            'where the angles are given in degrees', FutureWarning)
        self.euler_rotate(phi * 180 / pi, theta * 180 / pi, psi * 180 / pi,
                          center)

    def euler_rotate(self, phi=0.0, theta=0.0, psi=0.0, center=(0, 0, 0)):
        """Rotate atoms via Euler angles (in degrees).

        See e.g http://mathworld.wolfram.com/EulerAngles.html for explanation.

        Parameters:

        center :
            The point to rotate about. A sequence of length 3 with the
            coordinates, or 'COM' to select the center of mass, 'COP' to
            select center of positions or 'COU' to select center of cell.
        phi :
            The 1st rotation angle around the z axis.
        theta :
            Rotation around the x axis.
        psi :
            2nd rotation around the z axis.

        """
        if isinstance(center, str):
            if center.lower() == 'com':
                center = self.get_center_of_mass()
            elif center.lower() == 'cop':
                center = self.get_positions().mean(axis=0)
            elif center.lower() == 'cou':
                center = self.get_cell().sum(axis=0) / 2
            else:
                raise ValueError('Cannot interpret center')
        else:
            center = np.array(center)

        phi *= pi / 180
        theta *= pi / 180
        psi *= pi / 180

        # First move the molecule to the origin In contrast to MATLAB,
        # numpy broadcasts the smaller array to the larger row-wise,
        # so there is no need to play with the Kronecker product.
        rcoords = self.positions - center
        # First Euler rotation about z in matrix form
        D = np.array(((cos(phi), sin(phi), 0.),
                      (-sin(phi), cos(phi), 0.),
                      (0., 0., 1.)))
        # Second Euler rotation about x:
        C = np.array(((1., 0., 0.),
                      (0., cos(theta), sin(theta)),
                      (0., -sin(theta), cos(theta))))
        # Third Euler rotation, 2nd rotation about z:
        B = np.array(((cos(psi), sin(psi), 0.),
                      (-sin(psi), cos(psi), 0.),
                      (0., 0., 1.)))
        # Total Euler rotation
        A = np.dot(B, np.dot(C, D))
        # Do the rotation
        rcoords = np.dot(A, np.transpose(rcoords))
        # Move back to the rotation point
        self.positions = np.transpose(rcoords) + center

    def get_dihedral(self, a1, a2=None, a3=None, a4=None, mic=False):
        """Calculate dihedral angle.

        Calculate dihedral angle (in degrees) between the vectors a1->a2
        and a3->a4.

        Use mic=True to use the Minimum Image Convention and calculate the
        angle across periodic boundaries.
        """

        if a2 is None:
            # Old way - use radians
            warnings.warn(
                'Please_koopmans use new API (which will return the angle in degrees): '
                'atoms_obj.get_dihedral(a1,a2,a3,a4)*pi/180 instead of '
                'atoms_obj.get_dihedral([a1,a2,a3,a4])', FutureWarning)
            assert a3 is None and a4 is None
            a1, a2, a3, a4 = a1
            f = pi / 180
        else:
            f = 1

        if any(a is None for a in [a2, a3, a4]):
            raise ValueError('a2, a3 and a4 must not be None')

        # vector 1->2, 2->3, 3->4 and their normalized cross products:
        a = self.positions[a2] - self.positions[a1]
        b = self.positions[a3] - self.positions[a2]
        c = self.positions[a4] - self.positions[a3]
        if mic:
            a, b, c = find_mic([a, b, c], self.cell, self.pbc)[0]
        bxa = np.cross(b, a)
        cxb = np.cross(c, b)
        bxanorm = np.linalg.norm(bxa)
        cxbnorm = np.linalg.norm(cxb)
        if bxanorm == 0 or cxbnorm == 0:
            raise ZeroDivisionError('Undefined dihedral angle')
        bxa /= bxanorm
        cxb /= cxbnorm
        angle = np.vdot(bxa, cxb)
        # check for numerical trouble due to finite precision:
        if angle < -1:
            angle = -1
        if angle > 1:
            angle = 1
        angle = np.arccos(angle) * 180 / pi
        if np.vdot(bxa, c) > 0:
            angle = 360 - angle
        return angle * f

    def _masked_rotate(self, center, axis, diff, mask):
        # do rotation of subgroup by copying it to temporary atoms object
        # and then rotating that
        #
        # recursive object definition might not be the most elegant thing,
        # more generally useful might be a rotation function with a mask?
        group = self.__class__()
        for i in range(len(self)):
            if mask[i]:
                group += self[i]
        group.translate(-center)
        group.rotate(diff * 180 / pi, axis)
        group.translate(center)
        # set positions in original atoms object
        j = 0
        for i in range(len(self)):
            if mask[i]:
                self.positions[i] = group[j].position
                j += 1

    def set_dihedral(self, a1, a2=None, a3=None, a4=None, angle=None,
                     mask=None, indices=None):
        """Set the dihedral angle (degrees) between vectors a1->a2 and
        a3->a4 by changing the atom indexed by a4.

        If mask is not None, all the atoms described in mask
        (read: the entire subgroup) are moved. Alternatively to the mask,
        the indices of the atoms to be rotated can be supplied. If both
        *mask* and *indices* are given, *indices* overwrites *mask*.

        **Important**: If *mask* or *indices* is given and does not contain
        *a4*, *a4* will NOT be moved. In most case_koopmanss you therefore want
        to include *a4* in *mask*/*indices*.

        Example: the following defines a very crude
        ethane-like molecule and twists one half of it by 30 degrees.

        >>> atoms = Atoms('HHCCHH', [[-1, 1, 0], [-1, -1, 0], [0, 0, 0],
        ...                          [1, 0, 0], [2, 1, 0], [2, -1, 0]])
        >>> atoms.set_dihedral(1, 2, 3, 4, 210, mask=[0, 0, 0, 1, 1, 1])
        """

        if isinstance(a1, int):
            if any(a is None for a in [a2, a3, a4, angle]):
                raise ValueError('a2, a3, a4, and angle must not be None')
            angle *= pi / 180
        else:
            warnings.warn(
                'Please_koopmans use new API: '
                'atoms_obj.set_dihedral(a1,a2,a3,a4,angle) '
                'where angle is given in degrees', FutureWarning)
            if angle is None:
                angle = a2
                if mask is None:
                    mask = a3
                    if indices is None:
                        indices = a4
            else:
                assert a2 is None and a3 is None and a4 is None
            a1, a2, a3, a4 = a1

        # if not provided, set mask to the last atom in the
        # dihedral description
        if mask is None and indices is None:
            mask = np.zeros(len(self))
            mask[a4] = 1
        elif indices is not None:
            mask = [index in indices for index in range(len(self))]

        # compute necessary in dihedral change, from current value
        current = self.get_dihedral(a1, a2, a3, a4) * pi / 180
        diff = angle - current
        axis = self.positions[a3] - self.positions[a2]
        center = self.positions[a3]
        self._masked_rotate(center, axis, diff, mask)

    def rotate_dihedral(self, a1, a2=None, a3=None, a4=None,
                        angle=None, mask=None, indices=None):
        """Rotate dihedral angle.

        Same usage as in :meth:`ase_koopmans.Atoms.set_dihedral`: Rotate a group by a
        predefined dihedral angle, starting from its current configuration.
        """
        if not isinstance(a1, int):
            warnings.warn(
                'Please_koopmans use new API: '
                'atoms_obj.rotate_dihedral(a1,a2,a3,a4,angle) '
                'where angle is given in degrees', FutureWarning)
            if angle is None:
                angle = a2
                if mask is None and indices is None:
                    mask = a3
            else:
                assert a2 is None and a3 is None and a4 is None
            a1, a2, a3, a4 = a1

        if any(a is None for a in [a2, a3, a4, angle]):
            raise ValueError('a2, a3, a4, and angle must not be None')

        start = self.get_dihedral(a1, a2, a3, a4)
        self.set_dihedral(a1, a2, a3, a4, angle + start, mask, indices)

    def get_angle(self, a1, a2, a3, mic=False):
        """Get angle formed by three atoms.

        calculate angle in degrees between the vectors a2->a1 and
        a2->a3.

        Use mic=True to use the Minimum Image Convention and calculate the
        angle across periodic boundaries.
        """

        return self.get_angles([[a1, a2, a3]], mic=mic)[0]

    def get_angles(self, indices, mic=False):
        """Get angle formed by three atoms for multiple groupings.

        calculate angle in degrees between vectors between atoms a2->a1
        and a2->a3, where a1, a2, and a3 are in each row of indices.

        Use mic=True to use the Minimum Image Convention and calculate
        the angle across periodic boundaries.
        """

        indices = np.array(indices)

        a1s = self.positions[indices[:, 0]]
        a2s = self.positions[indices[:, 1]]
        a3s = self.positions[indices[:, 2]]

        v12 = a1s - a2s
        v32 = a3s - a2s

        cell = None
        pbc = None

        if mic:
            cell = self.cell
            pbc = self.pbc

        return get_angles(v12, v32, cell=cell, pbc=pbc)

    def set_angle(self, a1, a2=None, a3=None, angle=None, mask=None,
                  indices=None, add=False):
        """Set angle (in degrees) formed by three atoms.

        Sets the angle between vectors *a2*->*a1* and *a2*->*a3*.

        If *add* is `True`, the angle will be changed by the value given.

        Same usage as in :meth:`ase_koopmans.Atoms.set_dihedral`.
        If *mask* and *indices*
        are given, *indices* overwrites *mask*. If *mask* and *indices*
        are not set, only *a3* is moved."""

        if not isinstance(a1, int):
            # old API (uses radians)
            warnings.warn(
                'Please_koopmans use new API: '
                'atoms_obj.set_angle(a1,a2,a3,angle) '
                'where angle is given in degrees', FutureWarning)
            if angle is None:
                angle = a2
                if mask is None:
                    mask = a3
                a1, a2, a3 = a1
            else:
                assert a2 is None and a3 is None
            angle *= 180 / pi

        if any(a is None for a in [a2, a3, angle]):
            raise ValueError('a2, a3, and angle must not be None')

        # If not provided, set mask to the last atom in the angle description
        if mask is None and indices is None:
            mask = np.zeros(len(self))
            mask[a3] = 1
        elif indices is not None:
            mask = [index in indices for index in range(len(self))]

        if add:
            diff = angle
        else:
            # Compute necessary in angle change, from current value
            diff = angle - self.get_angle(a1, a2, a3)

        diff *= pi / 180
        # Do rotation of subgroup by copying it to temporary atoms object and
        # then rotating that
        v10 = self.positions[a1] - self.positions[a2]
        v12 = self.positions[a3] - self.positions[a2]
        v10 /= np.linalg.norm(v10)
        v12 /= np.linalg.norm(v12)
        axis = np.cross(v10, v12)
        center = self.positions[a2]
        self._masked_rotate(center, axis, diff, mask)

    def rattle(self, stdev=0.001, seed=None, rng=None):
        """Randomly displace atoms.

        This method adds random displacements to the atomic positions,
        taking a possible constraint into account.  The random numbers are
        drawn from a normal distribution of standard deviation stdev.

        For a parallel calculation, it is important to use the same
        seed on all processors!  """

        if seed is not None and rng is not None:
            raise ValueError('Please_koopmans do not provide both seed and rng.')

        if rng is None:
            if seed is None:
                seed = 42
            rng = np.random.RandomState(seed)
        positions = self.arrays['positions']
        self.set_positions(positions +
                           rng.normal(scale=stdev, size=positions.shape))

    def get_distance(self, a0, a1, mic=False, vector=False):
        """Return distance between two atoms.

        Use mic=True to use the Minimum Image Convention.
        vector=True gives the distance vector (from a0 to a1).
        """

        return self.get_distances(a0, [a1], mic=mic, vector=vector)[0]

    def get_distances(self, a, indices, mic=False, vector=False):
        """Return distances of atom No.i with a list of atoms.

        Use mic=True to use the Minimum Image Convention.
        vector=True gives the distance vector (from a to self[indices]).
        """

        R = self.arrays['positions']
        p1 = [R[a]]
        p2 = R[indices]

        cell = None
        pbc = None

        if mic:
            cell = self.cell
            pbc = self.pbc

        D, D_len = get_distances(p1, p2, cell=cell, pbc=pbc)

        if vector:
            D.shape = (-1, 3)
            return D
        else:
            D_len.shape = (-1,)
            return D_len

    def get_all_distances(self, mic=False, vector=False):
        """Return distances of all of the atoms with all of the atoms.

        Use mic=True to use the Minimum Image Convention.
        """
        R = self.arrays['positions']

        cell = None
        pbc = None

        if mic:
            cell = self.cell
            pbc = self.pbc

        D, D_len = get_distances(R, cell=cell, pbc=pbc)

        if vector:
            return D
        else:
            return D_len

    def set_distance(self, a0, a1, distance, fix=0.5, mic=False,
                     mask=None, indices=None, add=False, factor=False):
        """Set the distance between two atoms.

        Set the distance between atoms *a0* and *a1* to *distance*.
        By default, the center of the two atoms will be fixed.  Use
        *fix=0* to fix the first atom, *fix=1* to fix the second
        atom and *fix=0.5* (default) to fix the center of the bond.

        If *mask* or *indices* are set (*mask* overwrites *indices*),
        only the atoms defined there are moved
        (see :meth:`ase_koopmans.Atoms.set_dihedral`).

        When *add* is true, the distance is changed by the value given.
        In combination
        with *factor* True, the value given is a factor scaling the distance.

        It is assumed that the atoms in *mask*/*indices* move together
        with *a1*. If *fix=1*, only *a0* will therefore be moved."""

        if a0 % len(self) == a1 % len(self):
            raise ValueError('a0 and a1 must not be the same')

        if add:
            oldDist = self.get_distance(a0, a1, mic=mic)
            if factor:
                newDist = oldDist * distance
            else:
                newDist = oldDist + distance
            self.set_distance(a0, a1, newDist, fix=fix, mic=mic,
                              mask=mask, indices=indices, add=False,
                              factor=False)
            return

        R = self.arrays['positions']
        D = np.array([R[a1] - R[a0]])

        if mic:
            D, D_len = find_mic(D, self.cell, self.pbc)
        else:
            D_len = np.array([np.sqrt((D**2).sum())])
        x = 1.0 - distance / D_len[0]

        if mask is None and indices is None:
            indices = [a0, a1]
        elif mask:
            indices = [i for i in range(len(self)) if mask[i]]

        for i in indices:
            if i == a0:
                R[a0] += (x * fix) * D[0]
            else:
                R[i] -= (x * (1.0 - fix)) * D[0]

    def get_scaled_positions(self, wrap=True):
        """Get positions relative to unit cell.

        If wrap is True, atoms outside the unit cell will be wrapped into
        the cell in those directions with periodic boundary conditions
        so that the scaled coordinates are between zero and one."""

        fractional = self.cell.scaled_positions(self.positions)

        if wrap:
            for i, periodic in enumerate(self.pbc):
                if periodic:
                    # Yes, we need to do it twice.
                    # See the scaled_positions.py test.
                    fractional[:, i] %= 1.0
                    fractional[:, i] %= 1.0

        return fractional

    def set_scaled_positions(self, scaled):
        """Set positions relative to unit cell."""
        self.positions[:] = self.cell.cartesian_positions(scaled)

    def wrap(self, **wrap_kw):
        """Wrap positions to unit cell.

        Parameters:

        wrap_kw: (keyword=value) pairs
            optional keywords `pbc`, `center`, `pretty_translation`, `eps`,
            see :func:`ase_koopmans.geometry.wrap_positions`
        """

        if 'pbc' not in wrap_kw:
            wrap_kw['pbc'] = self.pbc

        self.positions[:] = self.get_positions(wrap=True, **wrap_kw)

    def get_temperature(self):
        """Get the temperature in Kelvin."""
        dof = len(self) * 3
        for constraint in self._constraints:
            dof -= constraint.removed_dof
        ekin = self.get_kinetic_energy()
        return 2 * ekin / (dof * units.kB)

    def __eq__(self, other):
        """Check for identity of two atoms objects.

        Identity means: same positions, atomic numbers, unit cell and
        periodic boundary conditions."""
        if not isinstance(other, Atoms):
            return False
        a = self.arrays
        b = other.arrays
        return (len(self) == len(other) and
                (a['positions'] == b['positions']).all() and
                (a['numbers'] == b['numbers']).all() and
                (self.cell == other.cell).all() and
                (self.pbc == other.pbc).all())

    def __ne__(self, other):
        """Check if two atoms objects are not equal.

        Any differences in positions, atomic numbers, unit cell or
        periodic boundary condtions make atoms objects not equal.
        """
        eq = self.__eq__(other)
        if eq is NotImplemented:
            return eq
        else:
            return not eq

    def get_volume(self):
        """Get volume of unit cell."""
        if self.cell.rank != 3:
            raise ValueError(
                'You have {0} lattice vectors: volume not defined'
                .format(self.cell.rank))
        return self.cell.volume

    def _get_positions(self):
        """Return reference to positions-array for in-place manipulations."""
        return self.arrays['positions']

    def _set_positions(self, pos):
        """Set positions directly, bypassing constraints."""
        self.arrays['positions'][:] = pos

    positions = property(_get_positions, _set_positions,
                         doc='Attribute for direct ' +
                         'manipulation of the positions.')

    @property
    def adsorbate_info(self):
        """Return the adsorbate information set by one of the surface
        builder functions. This function is only supplied in order to give
        a warning if this attribute (atoms.adsorbate_info) is asked for.
        The dictionary with adsorbate information has been moved to the
        info dictionary, i.e. atoms.info['adsorbate_info']."""
        warnings.warn("The adsorbate_info dictionary has been moved" +
                      " inside the info dictionary, i.e. atoms." +
                      "info['adsorbate_info']", FutureWarning)
        return self.info['adsorbate_info']

    @adsorbate_info.setter
    def adsorbate_info(self, dct):
        warnings.warn("The adsorbate_info dictionary has been moved" +
                      " inside the info dictionary, i.e. atoms." +
                      "info['adsorbate_info']", FutureWarning)
        self.info['adsorbate_info'] = dct

    def _get_atomic_numbers(self):
        """Return reference to atomic numbers for in-place
        manipulations."""
        return self.arrays['numbers']

    numbers = property(_get_atomic_numbers, set_atomic_numbers,
                       doc='Attribute for direct ' +
                       'manipulation of the atomic numbers.')

    def _get_cell(self):
        """Return reference to unit cell for in-place manipulations."""
        return self._cellobj

    cell = property(_get_cell, set_cell, doc='Attribute for direct ' +
                    'manipulation of the unit :class:`ase_koopmans.cell.Cell`.')

    def _get_pbc(self):
        """Return reference to pbc-flags for in-place manipulations."""
        return self._pbc

    pbc = property(_get_pbc, set_pbc,
                   doc='Attribute for direct manipulation ' +
                   'of the periodic boundary condition flags.')

    def write(self, filename, format=None, **kwargs):
        """Write atoms object to a file.

        see ase_koopmans.io.write for formats.
        kwargs are passed to ase_koopmans.io.write.
        """
        from ase_koopmans.io import write
        write(filename, self, format, **kwargs)

    def iterimages(self):
        yield self

    def edit(self):
        """Modify atoms interactively through ASE's GUI viewer.

        Conflicts leading to undesirable behaviour might arise
        when matplotlib has been pre-imported with certain
        incompatible backends and while trying to use the
        plot feature inside the interactive GUI. To circumvent,
        please_koopmans set matplotlib.use('gtk') before calling this
        method.
        """
        from ase_koopmans.gui.gui import GUI
        from ase_koopmans.gui.images import Images
        images = Images([self])
        gui = GUI(images)
        gui.run()


def string2vector(v):
    if isinstance(v, str):
        if v[0] == '-':
            return -string2vector(v[1:])
        w = np.zeros(3)
        w['xyz'.index(v)] = 1.0
        return w
    return np.array(v, float)


def default(data, dflt):
    """Helper function for setting default values."""
    if data is None:
        return None
    elif isinstance(data, (list, tuple)):
        newdata = []
        allnone = True
        for x in data:
            if x is None:
                newdata.append(dflt)
            else:
                newdata.append(x)
                allnone = False
        if allnone:
            return None
        return newdata
    else:
        return data
