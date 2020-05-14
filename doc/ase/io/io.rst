.. module:: ase.io
   :synopsis: File input-output module

=====================
File input and output
=====================

.. seealso::

    * :mod:`ase.io.trajectory`

.. toctree::
    :hidden:

    formatoptions
    trajectory
    ulm
    opls


The :mod:`ase.io` module has three basic functions: :func:`read`,
:func:`iread` and :func:`write`. The methods are described here:

.. autofunction:: read
.. autofunction:: iread
.. autofunction:: write

Use ``ase info --formats`` to see a list of formats.  This information
is programmatically accessible as ``ase.io.formats.ioformats``, a
dictionary which maps format names to :class:`ase.io.formats.IOFormat`
objects.

These are the file-formats that are recognized (formats with a ``+`` support
multiple configurations):

.. csv-table::
    :file: io.csv
    :header-rows: 1

.. note::

    Even though that ASE does a good job reading the above listed
    formats, it may not read some unusual features or strangely
    formatted files.

    For the CIF format, STAR extensions as save frames, global blocks,
    nested loops and multi-data values are not supported.

.. note::

    ASE read and write functions are automatically parallelized if a
    suitable MPI library is found. This requires to call read and write
    with same input on all cores. For more information, see
    :mod:`ase.parallel`.

.. note::

    ASE can read and write directly to compressed files. Simply add ``.gz``,
    ``.bz2`` or ``.xz`` to your filename.

The :func:`read` function is only designed to retrieve the atomic configuration
from a file, but for the CUBE format you can import the function:

.. function:: read_cube_data


which will return a ``(data, atoms)`` tuple::

  from ase.io.cube import read_cube_data
  data, atoms = read_cube_data('abc.cube')


Examples
========

>>> from ase import Atoms
>>> from ase.build import fcc111, add_adsorbate, bulk
>>> from ase.io import read, write
>>> adsorbate = Atoms('CO')
>>> adsorbate[1].z = 1.1
>>> a = 3.61
>>> slab = fcc111('Cu', (2, 2, 3), a=a, vacuum=7.0)
>>> add_adsorbate(slab, adsorbate, 1.8, 'ontop')

Write PNG image

>>> write('slab.png', slab * (3, 3, 1), rotation='10z,-80x')

.. image:: io1.png

Write animation with 500 ms duration per frame

>>> write('movie.gif', [bulk(s) for s in ['Cu', 'Ag', 'Au']], interval=500)


Write POVRAY file

>>> write('slab.pov', slab * (3, 3, 1), rotation='10z,-80x')

This will write both a ``slab.pov`` and a ``slab.ini`` file.  Convert
to PNG with the command ``povray slab.ini`` or use the
``run_povray=True`` option:

.. image:: io2.png

Here is an example using ``bbox``

>>> d = a / 2**0.5
>>> write('slab.pov', slab * (2, 2, 1),
...       bbox=(d, 0, 3 * d, d * 3**0.5))

.. image:: io3.png

This is an axample of display bond order for molecule

.. literalinclude:: save_C2H4.py

.. image:: C2H4.png

Note that in general the XYZ-format does not contain information about the unit cell, however, ASE uses the extended XYZ-format which stores the unitcell:

>>> from ase.io import read, write
>>> write('slab.xyz', slab)
>>> a = read('slab.xyz')
>>> cell = a.get_cell()
>>> cell.round(3)
array([[  5.105,   0.   ,   0.   ],
       [  2.553,   4.421,   0.   ],
       [  0.   ,   0.   ,  18.168]])
>>> a.get_pbc()
array([ True,  True, False], dtype=bool)

Another way to include the unit cell is to write the cell vectors at the end of the file as ``VEC<N> <x> <y> <z>`` (used for example in the ADF software).

>>> write('slab.xyz', vec_cell=True)

Use ASE's native format for writing all information:

>>> write('slab.traj', slab)
>>> b = read('slab.traj')
>>> b.cell.round(3)
array([[  5.105,   0.   ,   0.   ],
       [  2.553,   4.421,   0.   ],
       [  0.   ,   0.   ,  18.168]])
>>> b.pbc
array([ True,  True, False], dtype=bool)

A script showing all of the povray parameters, and generating the image below,
can be found here: :download:`save_pov.py`

.. image:: NaCl_C6H6.png

An other example showing how to change colors and textures in pov can
be found here: :download:`../../tutorials/saving_graphics.py`.

Adding a new file-format to ASE
===============================

Try to model the read/write functions after the *xyz* format as implemented
in :git:`ase/io/xyz.py` and also read, understand and update
:git:`ase/io/formats.py`.
