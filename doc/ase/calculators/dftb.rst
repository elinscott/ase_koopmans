.. module:: ase.calculators.dftb

=========
DftbPlus
=========

Introduction
============

DftbPlus_ is a density-functional based tight-binding code using
atom-centered orbitals. This interface makes it possible to use DftbPlus_
as a calculator in ASE. You need Slater-Koster files for the combination
of atom types of your system. These can be obtained at dftb.org_.

.. _DftbPlus: https://www.dftbplus.org/
.. _dftb.org: http://www.dftb.org/



Environment variables
=====================

The default command that ASE will use to start DftbPlus is
``dftb+ > PREFIX.out``. You can change this command by setting the
:envvar:`ASE_DFTB_COMMAND` environment variable, e.g.::

  $ export ASE_DFTB_COMMAND="/path/to/dftb+ > PREFIX.out"

For compatibility, also the old :envvar:`DFTB_COMMAND` variable can
be used, and the resulting command will be ``$DFTB_COMMAND > PREFIX.out``.
Before each DftbPlus calculation, also make sure that the
:envvar:`DFTB_PREFIX` variable points to the directory where
the Slater-Koster files are kept, e.g.::

  $ export DFTB_PREFIX=/path/to/mio-0-1/


DftbPlus Calculator (a FileIOCalculator)
========================================
The file 'geo_end.gen' contains the input and output geometry
and it will be updated during the DftbPlus calculations.

All keywords to the DftbPlus calculator can be set by ASE.


Parameters
==========
        restart: str (default None)
            If restart == None
            it is assumed that a new input file 'dftb_hsd.in'
            will be written by ase using default keywords
            and the ones given by the user.

            If restart != None
            it is assumed that keywords are in file 'restart'
        ignore_bad_restart_file: bool (default False)
            Ignore broken or missing restart file.  By default, an error
            is raised if restart!=None and the restart file is missing
            or broken.
        label: str (default 'dftb')
            Prefix used for the main output file (<label>.out).
        atoms: Atoms object (default None)
            Optional Atoms object to which the calculator will be attached.
            When restarting, atoms will get its positions and unit cell updated
            from the restart file.
        kpts: (default None)
            Brillouin zone sampling:

            * ``(1,1,1)`` or ``None``: Gamma-point only
            * ``(n1,n2,n3)``: Monkhorst-Pack grid
            * ``dict``: Interpreted as a path in the Brillouin zone if it
              contains the 'path_' keyword. Otherwise it is converted into a
              Monkhorst-Pack grid using ase.calculators.calculator.kpts2sizeandoffsets
            * ``[(k11,k12,k13),(k21,k22,k23),...]``: Explicit (Nkpts x 3) array of k-points
              in units of the reciprocal lattice vectors (each with equal weight)

.. _path: https://wiki.fysik.dtu.dk/ase/ase/dft/kpoints.html#ase.dft.kpoints.bandpath


Example: Geometry Optimization by ASE
=====================================

.. literalinclude:: dftb_ex1_relax.py

Example: Geometry Optimization by DftbPlus
==========================================

.. literalinclude:: dftb_ex2_relaxbyDFTB.py

Example: NVE md followed by NVT md (both by DftbPlus)
=====================================================

This is unphysical because of at least two reasons:

- oxygen does not have spin here
- the Berendsen coupling is too strong (0.01 here should be 0.0001)

.. literalinclude:: dftb_ex3_make_2h2o.py



