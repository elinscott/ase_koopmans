Geometry tools
==============

.. automodule:: ase.geometry
    :members:

Analysis tools
--------------
.. currentmodule:: ase.geometry.analysis

Provides the class :class:`Analysis` for structural analysis of any :class:`~ase.Atoms` object or list thereof (trajectories).

Example:

>>> import numpy as np
>>> from ase.build import molecule
>>> from ase.geometry.analysis import Analysis
>>> mol = molecule('C60')
>>> ana = Analysis(mol)
>>> CCBonds = ana.get_bonds('C', 'C', unique=True)
>>> CCCAngles = ana.get_angles('C', 'C', 'C', unique=True)
>>> print("There are {} C-C bonds in C60.".format(len(CCBonds[0])))
>>> print("There are {} C-C-C angles in C60.".format(len(CCCAngles[0])))
>>> CCBondValues = ana.get_values(CCBonds)
>>> CCCAngleValues = ana.get_values(CCCAngles)
>>> print("The average C-C bond length is {}.".format(np.average(CCBondValues)))
>>> print("The average C-C-C angle is {}.".format(np.average(CCCAngleValues)))

The :class:`Analysis` class provides a getter and setter for the images.
This allows you to use the same neighbourlist for different images, e.g. to analyze two MD simulations at different termperatures but constant bonding patterns.
Using this approach saves the time to recalculate all bonds, angles and dihedrals and therefore speeds up your analysis.

Using the :func:`Analysis.clear_cache()` function allows you to clear the calculated matrices/lists to reduce your memory usage.

The entire class can be used with few commands:

* To retrieve tuples of bonds/angles/dihedrals (they are calculated the first time they are accessed) use ``instance.all_xxx`` where *xxx* is one of bonds/angles/dihedrals.
* If you only want those one-way (meaning e.g. not bonds i-j and j-i but just i-j) use ``instance.unique_xxx``.
* To get selected bonds/angles/dihedrals use ``instance.get_xxx(A,B,...)``, see the API section for details on which arguments you can pass.
* To get the actual value of a bond/angle/dihedral use ``instance.get_xxx_value(tuple)``.
* To get a lot of bond/angle/dihedral values at once use :func:`Analysis.get_values()`.
* There is also a wrapper to get radial distribution functions :func:`Analysis.get_rdf()`.

The main difference between properties (getters) and functions here is, that getters provide data that is cached.
This means that getting information from ``Analysis.all_bonds`` more than once is instantaneous, since the information is cached in ``Analysis._cache``.
If you call any ``Analysis.get_xxx()`` the information is calculated from the cached data, meaning each call will take the same amount of time.


**API:**

.. currentmodule:: ase.geometry.analysis.Analysis
.. autoclass:: ase.geometry.analysis.Analysis
    :members:
