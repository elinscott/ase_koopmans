.. module:: ase.dft.kpoints
   :synopsis: Brillouin zone sampling

=======================
Brillouin zone sampling
=======================

The **k**-points are always given relative to the basis vectors of the
reciprocal unit cell.


Monkhorst-Pack
--------------

.. autofunction:: monkhorst_pack

The k-points are given as [MonkhorstPack]_:

.. math::

    \sum_{i=1,2,3} \frac{2n_i -N_i - 1}{2N_i} \mathbf{b}_i,

where `n_i=1,2,...,N_i`, ``size`` = `(N_1, N_2, N_3)` and the
`\mathbf{b}_i`'s are reciprocal lattice vectors.

.. autofunction:: get_monkhorst_pack_size_and_offset

Example:

>>> from ase.dft.kpoints import *
>>> monkhorst_pack((4, 1, 1))
array([[-0.375,  0.   ,  0.   ],
       [-0.125,  0.   ,  0.   ],
       [ 0.125,  0.   ,  0.   ],
       [ 0.375,  0.   ,  0.   ]])
>>> get_monkhorst_pack_size_and_offset([[0, 0, 0]])
(array([1, 1, 1]), array([ 0.,  0.,  0.]))


.. [MonkhorstPack]
    Hendrik J. Monkhorst and James D. Pack:
    *Special points for Brillouin-zone integrations*,
    Phys. Rev. B 13, 5188–5192 (1976)


Special points in the Brillouin zone
------------------------------------

.. data:: special_points

The below table lists the special points from [Setyawan-Curtarolo]_.


.. toctree:: bztable

.. include:: bztable.rst

.. [Setyawan-Curtarolo]
    High-throughput electronic band structure calculations:
    Challenges and tools

    Wahyu Setyawan, Stefano Curtarolo

    Computational Materials Science,
    Volume 49, Issue 2, August 2010, Pages 299–312

    https://doi.org/10.1016/j.commatsci.2010.05.010

You can find the special points in the Brillouin zone:

>>> from ase.build import bulk
>>> si = bulk('Si', 'diamond', a=5.459)
>>> lat = si.cell.get_bravais_lattice()
>>> print(list(lat.get_special_points()))
['G', 'K', 'L', 'U', 'W', 'X']
>>> path = si.cell.bandpath('GXW', npoints=100)
>>> print(path.kpts.shape)
(100, 3)

.. autofunction:: get_special_points
.. autofunction:: bandpath
.. autofunction:: parse_path_string
.. autofunction:: labels_from_kpts


Band path
---------

The :class:`~ase.dft.kpoints.BandPath` class stores all the relevant
band path information in a single object.
It is typically created by helper functions such as
:meth:`ase.cell.Cell.bandpath` or :meth:`ase.lattice.BravaisLattice.bandpath`.

.. autoclass:: BandPath
               :members:

Band structure
--------------

.. autoclass:: ase.dft.band_structure.BandStructure
   :members:

Free electron example:

.. literalinclude:: bs.py

.. image:: cu.png


Interpolation
-------------

.. autofunction:: monkhorst_pack_interpolate


High symmetry paths
-------------------

.. data:: special_paths

The :mod:`ase.lattice` framework provides suggestions for high symmetry
paths in the BZ from the [Setyawan-Curtarolo]_ paper.

>>> from ase.lattice import BCC
>>> lat = BCC(3.5)
>>> lat.get_special_points()
{'G': array([0, 0, 0]), 'H': array([ 0.5, -0.5,  0.5]), 'P': array([0.25, 0.25, 0.25]), 'N': array([0. , 0. , 0.5])}
>>> lat.special_path
'GHNGPH,PN'

In case you want this information *without* providing the lattice parameter
(which is possible for those lattices where the points do not depend on the
lattice parameters), the data is also available as dictionaries:

>>> from ase.dft.kpoints(import special_paths, sc_special_points,
...                      parse_path_string)
>>> paths = sc_special_paths['bcc']
>>> paths
[['G', 'H', 'N', 'G', 'P', 'H'], ['P', 'N']]
>>> points = sc_special_points['bcc']
>>> points
{'H': [0.5, -0.5, 0.5], 'N': [0, 0, 0.5], 'P': [0.25, 0.25, 0.25],
 'G': [0, 0, 0]}
>>> kpts = [points[k] for k in paths[0]]  # G-H-N-G-P-H
>>> kpts
[[0, 0, 0], [0.5, -0.5, 0.5], [0, 0, 0.5], [0, 0, 0], [0.25, 0.25, 0.25], [0.5, -0.5, 0.5]]


Chadi-Cohen
-----------

Predefined sets of **k**-points:

.. data:: cc6_1x1
.. data:: cc12_2x3
.. data:: cc18_sq3xsq3
.. data:: cc18_1x1
.. data:: cc54_sq3xsq3
.. data:: cc54_1x1
.. data:: cc162_sq3xsq3
.. data:: cc162_1x1


Naming convention: ``cc18_sq3xsq3`` is 18 **k**-points for a
sq(3)xsq(3) cell.

Try this:

>>> import numpy as np
>>> import matplotlib.pyplot as plt
>>> from ase.dft.kpoints import cc162_1x1
>>> B = [(1, 0, 0), (-0.5, 3**0.5 / 2, 0), (0, 0, 1)]
>>> k = np.dot(cc162_1x1, B)
>>> plt.plot(k[:, 0], k[:, 1], 'o')  # doctest: +SKIP
[<matplotlib.lines.Line2D object at 0x9b61dcc>]
>>> plt.show()

.. image:: cc.png
