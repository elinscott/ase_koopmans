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

**API:**

.. currentmodule:: ase.geometry.analysis.Analysis
.. autoclass:: ase.geometry.analysis.Analysis
    :members:
