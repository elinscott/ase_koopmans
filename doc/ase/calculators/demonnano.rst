.. module:: ase.calculators.demonnano

==========
deMon-Nano
==========

Introduction
============

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

The file 'deMon.inp' contains the input geometry and parameters.
The file 'deMon.out' contains the results.

Parameters
==========

        label : str 
            relative path to the run directory
        atoms  : Atoms object
            the atoms object
        command  : str
            Command to run deMon. If not present, the environment variable ASE_DEMONNANO_COMMAND is used
        basis_path  : str 
            Relative path to the directory containing DFTB-SCC or DFTB-0 parameters
            If not present, the environment variable DEMONNANO_BASIS_PATH is used
        restart_path  : str 
            Relative path to the deMon restart dir
        title : str 
            Title in the deMon input file.
        forces : bool
            If True a force calculation is enforced
        print_out : str | list 
            Options for the printing in deMon
        input_arguments : dict 
            Explicitly given input arguments. The key is the input keyword
            and the value is either a str, a list of str (will be written on the same line as the keyword),
            or a list of lists of str (first list is written on the first line, the others on following lines.)

Example: Geometry Optimization with ASE
=======================================

.. literalinclude:: demonnano_ex1_relax.py

