.. _eos:

=======================
Equation of state (EOS)
=======================

First, do a bulk calculation for different lattice constants:

.. literalinclude:: eos1.py

This will write a trajectory file containing five configurations of
FCC silver for five different lattice constants.  Now, analyse the
result with the :class:`~ase.eos.EquationOfState` class and this
script:

.. literalinclude:: eos2.py

|eos|

A quicker way to do this analysis, is to use the :mod:`ase.gui` tool:

.. highlight:: bash

::

    $ ase gui Ag.traj

And then choose :menuselection:`Tools --> Bulk modulus`.

.. |eos| image:: Ag-eos.png
