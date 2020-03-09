.. _dimtutorial:

=======================
Dimensionality analysis
=======================

This is a example of analysis of the dimensionality of a structure using
the :func:`ase.geometry.dimensionality.analyze_dimensionality` function. This is
useful for finding low-dimensional materials, such as 1D chain-like
structures, 2D layered structures, or structures with multiple dimensionality
types, such as 1D+3D.

The example below creates a layered :mol:`MoS_2` structure and analyzes its
dimensionality.

.. literalinclude:: dimexample.py

Coloring the atoms by their tags shows the distinct bonded clusters, which in
this case are separate layers.

Each component in the material can be extracted, or "*isolated*",
using the :func:`ase.geometry.dimensionality.isolate_components` function as
the example below demonstrates.

.. literalinclude:: isolation_example.py

The method is described in the article:

  | P.M. Larsen, M. Pandey, M. Strange, and K. W. Jacobsen
  | `Definition of a scoring parameter to identify low-dimensional materials components`__
  | Phys. Rev. Materials 3 034003, 2019

__ https://doi.org/10.1103/PhysRevMaterials.3.034003

A preprint is available `here <https://arxiv.org/pdf/1808.02114.pdf>`_.

.. seealso::

    More examples here: `Dimensionality analysis of ICSD and COD databases
    <https://cmr.fysik.dtu.dk/lowdim/lowdim.html>`_.


.. autofunction:: ase.geometry.dimensionality.analyze_dimensionality

.. autofunction:: ase.geometry.dimensionality.isolate_components
