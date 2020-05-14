.. module:: ase.atoms
.. module:: ase

================
The Atoms object
================

The :class:`Atoms` object is a collection of atoms.  Here
is how to define a CO molecule::

  from ase import Atoms
  d = 1.1
  co = Atoms('CO', positions=[(0, 0, 0), (0, 0, d)])

Here, the first argument specifies the type of the atoms and we used
the ``positions`` keywords to specify their positions.  Other
possible keywords are: ``numbers``, ``tags``, ``momenta``, ``masses``,
``magmoms`` and ``charges``.

Here is how you could make an infinite gold wire with a bond length of
2.9 Å::

  from ase import Atoms
  d = 2.9
  L = 10.0
  wire = Atoms('Au',
               positions=[[0, L / 2, L / 2]],
               cell=[d, L, L],
               pbc=[1, 0, 0])

.. image:: Au-wire.png

Here, two more optional keyword arguments were used:

``cell``: Unit cell size
  This can be a sequence of three numbers for
  an orthorhombic unit cell or three by three numbers for a general
  unit cell (a sequence of three sequences of three numbers) or six numbers
  (three legths and three angles in degrees).  The default value is
  *[0,0,0]* which is the same as
  *[[0,0,0],[0,0,0],[0,0,0]]* or *[0,0,0,90,90,90]* meaning that none of the
  three lattice vectors are defined.

``pbc``: Periodic boundary conditions
  The default value is *False* - a value of *True* would give
  periodic boundary conditions along all three axes.  It is possible
  to give a sequence of three booleans to specify periodicity along
  specific axes.

You can also use the following methods to work with the unit cell and
the boundary conditions: :meth:`~Atoms.set_pbc`,
:meth:`~Atoms.set_cell`, :meth:`~Atoms.get_cell`,
and :meth:`~Atoms.get_pbc`.


Working with the array methods of Atoms objects
===============================================

Like with a single :class:`~ase.atom.Atom` the properties of a collection of atoms
can be accessed and changed with get- and set-methods. For example
the positions of the atoms can be addressed as

>>> from ase import Atoms
>>> atoms = Atoms('N3', [(0, 0, 0), (1, 0, 0), (0, 0, 1)])
>>> atoms.get_positions()
array([[ 0.,  0.,  0.],
       [ 1.,  0.,  0.],
       [ 0.,  0.,  1.]])
>>> atoms.set_positions([(2, 0, 0), (0, 2, 2), (2, 2, 0)])
>>> atoms.get_positions()
array([[ 2.,  0.,  0.],
       [ 0.,  2.,  2.],
       [ 2.,  2.,  0.]])

Here is the full list of the get/set methods operating on all the
atoms at once.  The get methods return an array of quantities, one for
each atom; the set methods take similar arrays.
E.g. :meth:`~Atoms.get_positions` return N * 3 numbers,
:meth:`~Atoms.get_atomic_numbers` return N integers.

*These methods return copies of the internal arrays.  It is thus safe
to modify the returned arrays.*

.. list-table::

  * - :meth:`~Atoms.get_atomic_numbers`
    - :meth:`~Atoms.set_atomic_numbers`
  * - :meth:`~Atoms.get_initial_charges`
    - :meth:`~Atoms.set_initial_charges`
  * - :meth:`~Atoms.get_charges`
    -
  * - :meth:`~Atoms.get_chemical_symbols`
    - :meth:`~Atoms.set_chemical_symbols`
  * - :meth:`~Atoms.get_initial_magnetic_moments`
    - :meth:`~Atoms.set_initial_magnetic_moments`
  * - :meth:`~Atoms.get_magnetic_moments`
    -
  * - :meth:`~Atoms.get_masses`
    - :meth:`~Atoms.set_masses`
  * - :meth:`~Atoms.get_momenta`
    - :meth:`~Atoms.set_momenta`
  * - :meth:`~Atoms.get_forces`
    -
  * - :meth:`~Atoms.get_positions`
    - :meth:`~Atoms.set_positions`
  * - :meth:`~Atoms.get_potential_energies`
    -
  * - :meth:`~Atoms.get_scaled_positions`
    - :meth:`~Atoms.set_scaled_positions`
  * - :meth:`~Atoms.get_stresses`
    -
  * - :meth:`~Atoms.get_tags`
    - :meth:`~Atoms.set_tags`
  * - :meth:`~Atoms.get_velocities`
    - :meth:`~Atoms.set_velocities`

There are also a number of get/set methods that operate on quantities
common to all the atoms or defined for the collection of atoms:

.. list-table::

  * - :meth:`~Atoms.get_calculator`
    - :meth:`~Atoms.set_calculator`
  * - :meth:`~Atoms.get_cell`
    - :meth:`~Atoms.set_cell`
  * - :meth:`~Atoms.get_cell_lengths_and_angles`
    -
  * - :meth:`~Atoms.get_center_of_mass`
    -
  * - :meth:`~Atoms.get_kinetic_energy`
    -
  * - :meth:`~Atoms.get_magnetic_moment`
    -
  * - :meth:`~Atoms.get_global_number_of_atoms`
    -
  * - :meth:`~Atoms.get_pbc`
    - :meth:`~Atoms.set_pbc`
  * - :meth:`~Atoms.get_potential_energy`
    -
  * - :meth:`~Atoms.get_stress`
    -
  * - :meth:`~Atoms.get_total_energy`
    -
  * - :meth:`~Atoms.get_volume`
    -


Unit cell and boundary conditions
=================================

The :class:`Atoms` object holds a unit cell.  The unit cell
is a :class:`~ase.cell.Cell` object which resembles a 3x3 matrix
when used with numpy, arithmetic operations, or indexing:

>>> atoms.cell
Cell([0.0, 0.0, 0.0], pbc=False)
>>> atoms.cell[:]
array([[0., 0., 0.],
       [0., 0., 0.],
       [0., 0., 0.]])

The cell can be defined or changed using the
:meth:`~Atoms.set_cell` method. Changing the unit cell
does per default not move the atoms:

>>> import numpy as np
>>> atoms.set_cell(2 * np.identity(3))
>>> atoms.get_cell()
Cell([2.0, 2.0, 2.0], pbc=False)
>>> atoms.set_positions([(2, 0, 0), (1, 1, 0), (2, 2, 0)])
>>> atoms.get_positions()
array([[ 2.,  0.,  0.],
       [ 1.,  1.,  0.],
       [ 2.,  2.,  0.]])

However if we set ``scale_atoms=True`` the atomic positions are scaled with
the unit cell:

>>> atoms.set_cell(np.identity(3), scale_atoms=True)
>>> atoms.get_positions()
array([[ 1. ,  0. ,  0. ],
       [ 0.5,  0.5,  0. ],
       [ 1. ,  1. ,  0. ]])

The :meth:`~Atoms.set_pbc` method specifies whether
periodic boundary conditions are to be used in the directions of the
three vectors of the unit cell.  A slab calculation with periodic
boundary conditions in *x* and *y* directions and free boundary
conditions in the *z* direction is obtained through

>>> atoms.set_pbc((True, True, False))

or

>>> atoms.pbc = (True, True, False)

.. _atoms_special_attributes:

Special attributes
==================

It is also possible to work directly with the attributes
:attr:`~Atoms.positions`, :attr:`~Atoms.numbers`,
:attr:`~Atoms.pbc` and :attr:`~Atoms.cell`.  Here
we change the position of the 2nd atom (which has count number 1
because Python starts counting at zero) and the type of the first
atom:

>>> atoms.positions *= 2
>>> atoms.positions[1] = (1, 1, 0)
>>> atoms.get_positions()
array([[ 2.,  0.,  0.],
       [ 1.,  1.,  0.],
       [ 2.,  2.,  0.]])
>>> atoms.positions
array([[ 2.,  0.,  0.],
       [ 1.,  1.,  0.],
       [ 2.,  2.,  0.]])
>>> atoms.numbers
array([7, 7, 7])
>>> atoms.numbers[0] = 13
>>> atoms.get_chemical_symbols()
['Al', 'N', 'N']

The atomic numbers can also be edited using the :attr:`~Atoms.symbols`
shortcut (see also :class:`ase.symbols.Symbols`):

>>> atoms.symbols
Symbols('AlN2')
>>> atoms.symbols[2] = 'Cu'
>>> atoms.symbols
Symbols('AlNCu')
>>> atoms.numbers
array([13,  7, 29])

Check for periodic boundary conditions:

>>> atoms.pbc  # equivalent to atoms.get_pbc()
array([ True,  True, False], dtype=bool)
>>> atoms.pbc.any()
True
>>> atoms.pbc[2] = 1
>>> atoms.pbc
array([ True,  True,  True], dtype=bool)

Hexagonal unit cell:

>>> atoms.cell = [2.5, 2.5, 15, 90, 90, 120]

Attributes that can be edited directly are:

* :meth:`~Atoms.numbers`
* :meth:`~Atoms.symbols`
* :meth:`~Atoms.positions`
* :meth:`~Atoms.cell`
* :meth:`~Atoms.pbc`
* :meth:`~Atoms.constraints`



Adding a calculator
===================

A calculator can be attached to the atoms with the purpose
of calculating energies and forces on the atoms. ASE works with many
different :mod:`ase.calculators`.

A calculator object *calc* is attached to the atoms like this:

>>> atoms.calc = calc

After the calculator has been appropriately setup the energy of the
atoms can be obtained through

>>> atoms.get_potential_energy()

The term "potential energy" here means for example the total energy of
a DFT calculation, which includes both kinetic, electrostatic, and
exchange-correlation energy for the electrons. The reason it is called
potential energy is that the atoms might also have a kinetic energy
(from the moving nuclei) and that is obtained with

>>> atoms.get_kinetic_energy()

In case of a DFT calculator, it is up to the user to check exactly what
the :meth:`~Atoms.get_potential_energy` method returns. For
example it may be the result of a calculation with a finite
temperature smearing of the occupation numbers extrapolated to zero
temperature.  More about this can be found for the different
:mod:`ase.calculators`.

The following methods can only be called if a calculator is present:

* :meth:`~Atoms.get_potential_energy`
* :meth:`~Atoms.get_potential_energies`
* :meth:`~Atoms.get_forces`
* :meth:`~Atoms.get_stress`
* :meth:`~Atoms.get_stresses`
* :meth:`~Atoms.get_total_energy`
* :meth:`~Atoms.get_magnetic_moments`
* :meth:`~Atoms.get_magnetic_moment`

Not all of these methods are supported by all calculators.


List-methods
============

.. list-table::

  * - method
    - example
  * - ``+``
    - ``wire2 = wire + co``
  * - ``+=``, :meth:`~Atoms.extend`
    - ``wire += co``

      ``wire.extend(co)``
  * - :meth:`~Atoms.append`
    - ``wire.append(Atom('H'))``
  * - ``*``
    - ``wire3 = wire * (3, 1, 1)``
  * - ``*=``, :meth:`~Atoms.repeat`
    - ``wire *= (3, 1, 1)``

      ``wire.repeat((3, 1, 1))``
  * - ``len``
    - ``len(co)``
  * - ``del``
    - ``del wire3[0]``

      ``del wire3[[1,3]]``
  * - :meth:`~Atoms.pop`
    - ``oxygen = wire2.pop()``


Note that the ``del`` method can be used with the more powerful numpy-style indexing, as in the second example above. This can be combined with python list comprehension in order to selectively delete atoms within an ASE Atoms object. For example, the below code creates an ethanol molecule and subsequently strips all the hydrogen atoms from it::

  from ase.build import molecule
  atoms = molecule('CH3CH2OH')
  del atoms[[atom.index for atom in atoms if atom.symbol=='H']]


Other methods
=============

* :meth:`~Atoms.center`
* :meth:`~Atoms.wrap`
* :meth:`~Atoms.translate`
* :meth:`~Atoms.rotate`
* :meth:`~Atoms.euler_rotate`
* :meth:`~Atoms.get_dihedral`
* :meth:`~Atoms.set_dihedral`
* :meth:`~Atoms.rotate_dihedral`
* :meth:`~Atoms.rattle`
* :meth:`~Atoms.set_constraint`
* :meth:`~Atoms.set_distance`
* :meth:`~Atoms.copy`
* :meth:`~Atoms.get_center_of_mass`
* :meth:`~Atoms.get_distance`
* :meth:`~Atoms.get_distances`
* :meth:`~Atoms.get_all_distances`
* :meth:`~Atoms.get_volume`
* :meth:`~Atoms.has`
* :meth:`~Atoms.edit`



List of all Methods
===================

.. autoclass:: Atoms
   :members:
