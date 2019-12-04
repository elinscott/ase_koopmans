.. module:: ase.calculators.demonnano

==========
deMon-Nano
==========

deMon-Nano_ is a density-functional based tight-binding (DFTB) code using atom centered orbitals. This
interface makes it possible to use deMon-Nano_ as a calculator in ASE.
You need Slater-Koster files for the combination of
atom types of your system. These can be obtained at dftb.org_.

.. _deMon-Nano: http://demon-nano.ups-tlse.fr/
.. _dftb.org: http://www.dftb.org/

Environment variables
=====================

Set environment variables in your configuration file (what is the directory
for the Slater-Koster files and what is the name of the executable):

- bash::

  $ DEMONNANO_BASIS_PATH="/path/to/basis/"  (an example)
  $ ASE_DEMONNANO_COMMAND="/path/to/bin/deMon.username.x (an example)

deMon-Nano Calculator (a FileIOCalculator)
==========================================

.. autoclass:: DemonNano

Example: Geometry Optimization with ASE
=======================================

.. literalinclude:: demonnano_ex1_relax.py

