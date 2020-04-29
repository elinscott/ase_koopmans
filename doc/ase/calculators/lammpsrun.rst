.. module:: ase.calculators.LAMMPSrun

=========
LAMMPSrun
=========

Introduction
============

LAMMPS_ (Large-scale Atomic/Molecular Massively Parallel Simulator) is a classical molecular dynamics code.

    "LAMMPS has potentials for soft materials (biomolecules, polymers) and solid-state materials (metals, semiconductors) and coarse-grained or mesoscopic systems. It can be used to model atoms or, more generically, as a parallel particle simulator at the atomic, meso, or continuum scale."


.. _LAMMPS: https://lammps.sandia.gov/

This is LAMMPSrun ASE implementation of the interface to LAMMPSRUN_.

Environment variables
=====================

The environment variable :envvar:`ASE_LAMMPSRUN_COMMAND` should contain
the path to the lammps binary, or more generally, a command line 
possibly also including an MPI-launcher command.
For example (in a Bourne-shell compatible environment):

.. highlight:: bash
 
::

  $ export ASE_LAMMPSRUN_COMMAND=/path/to/lmp_binary

.. highlight:: python

or possibly something similar to

.. highlight:: bash
 
::

  $ export ASE_LAMMPSRUN_COMMAND="/path/to/mpirun --np 4 lmp_binary"

.. highlight:: python



LAMMPS Calculator
================= 

The LAMMPSRUN calculator first appeared in ASE version 3.5.0.  At the
time of the release of ASE 3.17.0, the LAMMPS calculator is still a
thin wrapper containing basic features to enable the use of LAMMPS in
ASE (missing some feature might have been added in the source code
development tree or some more recent version of ASE).

.. class:: LAMMPS(..., keyword=value, ...)

Below follows a list with a selection of parameters

====================  ==========  ============== =============================
keyword               type        default value  description
====================  ==========  ============== =============================
``specorder``         ``list``    ``None``       List containing the atom
                                                 species;
``always_triclinic``  ``bool``    ``False``      LAMMPS treats orthorhombic and
                                                 tilted cells differently
                                                 (extra parameters in input
                                                 and output-files).
``keep_alive``        ``bool``    ``True``       Do not restart LAMMPS for
                                                 each query; perform only a
                                                 reset to a clean state.
``keep_tmp_files``    ``bool``    ``False``      Generated LAMMPS input and
                                                 output-files are not deleted.
``tmp_dir``           ``string``  ``None``       Location to write tmp-files.
``files``             ``list``    ``[]``         List of files needed by
                                                 LAMMPS. Typically a list of
                                                 potential files.
``verbose``           ``bool``    ``False``      Print additional debugging
                                                 output to STDOUT.
``write_velocities``  ``bool``    ``False``      Forward ASE velocities to
                                                 LAMMPS.
====================  ==========  ============== =============================



Example
=======

A simple example.

::

  from ase import Atoms, Atom
  from ase.calculators.lammpsrun import LAMMPS
  
  a = [6.5, 6.5, 7.7]
  d = 2.3608
  NaCl = Atoms([Atom('Na', [0, 0, 0]),
                Atom('Cl', [0, 0, d])],
               cell=a, pbc=True)
  
  calc = LAMMPS()
  NaCl.calc = calc
  
  print(NaCl.get_stress())

  
Setting up an OPLS calculation
==============================

There are some modules to facilitate the setup of an OPLS force field 
calculation, see :mod:`ase.io.opls`.


