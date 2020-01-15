.. module:: ase.cell

The Cell object
===============

The Cell object represents three lattice vectors forming a parallel epiped.

``atoms.cell`` is a Cell object.

Examples::

  >>> from ase.build import bulk
  >>> cell = bulk('Au').cell
  >>> cell
  Cell([[0.0, 2.04, 2.04], [2.04, 0.0, 2.04], [2.04, 2.04, 0.0]], pbc=True)

The cell behaves like a 3x3 array when used like one::

  >>> cell[:]
  array([[0.  , 2.04, 2.04],
         [2.04, 0.  , 2.04],
         [2.04, 2.04, 0.  ]])


Common functionality::

  >>> cell.lengths()
  array([2.88499567, 2.88499567, 2.88499567])
  >>> cell.angles()
  array([60., 60., 60.])
  >>> cell.volume
  16.979328000000002


.. autoclass:: Cell
    :members:
