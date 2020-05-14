.. module:: ase.calculators
   :synopsis: Energy, force and stress calculators.

.. _calculators:

===========
Calculators
===========

For ASE, a calculator is a black box that can take atomic numbers and
atomic positions from an :class:`~ase.Atoms` object and calculate the
energy and forces and sometimes also stresses.

In order to calculate forces and energies, you need to attach a
calculator object to your atoms object:

>>> atoms = read('molecule.xyz')
>>> e = atoms.get_potential_energy()  # doctest: IGNORE_EXCEPTION_DETAIL
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "/home/jjmo/ase/atoms/ase.py", line 399, in get_potential_energy
    raise RuntimeError('Atoms object has no calculator.')
RuntimeError: Atoms object has no calculator.
>>> from ase.calculators.abinit import Abinit
>>> calc = Abinit(...)
>>> atoms.calc = calc
>>> e = atoms.get_potential_energy()
>>> print(e)
-42.0

Here we attached
an instance of the :mod:`ase.calculators.abinit` class and then
we asked for the energy.


.. _supported calculators:

Supported calculators
=====================

The calculators can be divided in four groups:

1) Asap_, GPAW_, and Hotbit_ have their own native ASE interfaces.

2) ABINIT, AMBER, CP2K, CASTEP, deMon2k, DFTB+, ELK, EXCITING, FHI-aims, FLEUR, GAUSSIAN,
   Gromacs, LAMMPS, MOPAC, NWChem, Octopus, ONETEP, psi4, Q-Chem, Quantum ESPRESSO, SIESTA,
   TURBOMOLE and VASP, have Python wrappers in the ASE package, but the actual
   FORTRAN/C/C++ codes are not part of ASE.

3) Pure python implementations included in the ASE package: EMT, EAM,
   Lennard-Jones and Morse.

4) Calculators that wrap others, included in the ASE package:
   :class:`ase.calculators.checkpoint.CheckpointCalculator`,
   the :class:`ase.calculators.loggingcalc.LoggingCalculator`,
   the :class:`ase.calculators.socketio.SocketIOCalculator`,
   the :ref:`Grimme-D3 <grimme>` potential, and the qmmm calculators
   :class:`~ase.calculators.qmmm.EIQMMM`,  and :class:`~ase.calculators.qmmm.SimpleQMMM`.

========================================= ===========================================
name                                      description
========================================= ===========================================
Asap_                                     Highly efficient EMT code
GPAW_                                     Real-space/plane-wave/LCAO PAW code
Hotbit_                                   DFT based tight binding
:mod:`~ase.calculators.abinit`            Plane-wave pseudopotential code
:mod:`~ase.calculators.amber`             Classical molecular dynamics code
:mod:`~ase.calculators.castep`            Plane-wave pseudopotential code
:mod:`~ase.calculators.cp2k`              DFT and classical potentials
:mod:`~ase.calculators.demon`             Gaussian based DFT code
:mod:`~ase.calculators.demonnano`         DFT based tight binding code
:mod:`~ase.calculators.dftb`              DFT based tight binding
:mod:`~ase.calculators.dmol`              Atomic orbital DFT code
:mod:`~ase.calculators.eam`               Embedded Atom Method
elk                                       Full Potential LAPW code
:mod:`~ase.calculators.espresso`          Plane-wave pseudopotential code
:mod:`~ase.calculators.exciting`          Full Potential LAPW code
:mod:`~ase.calculators.aims`              Numeric atomic orbital, full potential code
:mod:`~ase.calculators.fleur`             Full Potential LAPW code
:mod:`~ase.calculators.gamess_us`         Gaussian based electronic structure code
:mod:`~ase.calculators.gaussian`          Gaussian based electronic structure code
:mod:`~ase.calculators.gromacs`           Classical molecular dynamics code
:mod:`~ase.calculators.gulp`              Interatomic potential code
:mod:`~ase.calculators.kim`               Classical MD with standardized models
:mod:`~ase.calculators.lammps`            Classical molecular dynamics code
:mod:`~ase.calculators.mixing`            Combination of multiple calculators
:mod:`~ase.calculators.mopac`             Semiempirical molecular orbital code
:mod:`~ase.calculators.nwchem`            Gaussian based electronic structure code
:mod:`~ase.calculators.octopus`           Real-space pseudopotential code
:mod:`~ase.calculators.onetep`            Linear-scaling pseudopotential code
:mod:`~ase.calculators.openmx`            LCAO pseudopotential code
:mod:`~ase.calculators.orca`              Gaussian based electronic structure code
:mod:`~ase.calculators.psi4`              Gaussian based electronic structure code
:mod:`~ase.calculators.qchem`             Gaussian based electronic structure code
:mod:`~ase.calculators.siesta`            LCAO pseudopotential code
:mod:`~ase.calculators.turbomole`         Fast atom orbital code
:mod:`~ase.calculators.vasp`              Plane-wave PAW code
:mod:`~ase.calculators.emt`               Effective Medium Theory calculator
lj                                        Lennard-Jones potential
morse                                     Morse potential
:mod:`~ase.calculators.checkpoint`        Checkpoint calculator
:mod:`~ase.calculators.socketio`          Socket-based interface to calculators
:mod:`~ase.calculators.loggingcalc`       Logging calculator
:mod:`~ase.calculators.dftd3`             DFT-D3 dispersion correction calculator
:class:`~ase.calculators.qmmm.EIQMMM`     Explicit Interaction QM/MM
:class:`~ase.calculators.qmmm.SimpleQMMM` Subtractive (ONIOM style) QM/MM
========================================= ===========================================

.. index:: D3, Grimme
.. _grimme:

.. note::

    A Fortran implemetation of the Grimme-D3 potential, that can be used as
    an add-on to any ASE calculator, can be found here:
    https://gitlab.com/ehermes/ased3/tree/master.

The calculators included in ASE are used like this:

>>> from ase.calculators.abc import ABC
>>> calc = ABC(...)

where ``abc`` is the module name and ``ABC`` is the class name.


.. _Asap: https://wiki.fysik.dtu.dk/asap
.. _GPAW: https://wiki.fysik.dtu.dk/gpaw
.. _Hotbit: https://github.com/pekkosk/hotbit

Calculator keywords
===================

Example for a hypothetical ABC calculator:

.. class:: ABC(restart=None, ignore_bad_restart_file=False, label=None,
               atoms=None, parameters=None, command='abc > PREFIX.abc',
               xc=None, kpts=[1, 1, 1], smearing=None,
               charge=0.0, nbands=None, **kwargs)

   Create ABC calculator

   restart: str
       Prefix for restart file.  May contain a directory.  Default
       is None: don't restart.
   ignore_bad_restart_file: bool
       Ignore broken or missing restart file.  By default, it is an
       error if the restart file is missing or broken.
   label: str
       Name used for all files.  May contain a directory.
   atoms: Atoms object
       Optional Atoms object to which the calculator will be
       attached.  When restarting, atoms will get its positions and
       unit-cell updated from file.
   command: str
       Command used to start calculation.  This will override any value
       in an :envvar:`ASE_ABC_COMMAND` environment variable.
   parameters: str
       Read parameters from file.
   xc: str
       XC-functional (``'LDA'``, ``'PBE'``, ...).
   kpts:
       Brillouin zone sampling:

       * ``(1,1,1)``: Gamma-point
       * ``(n1,n2,n3)``: Monkhorst-Pack grid
       * ``(n1,n2,n3,'gamma')``: Shifted Monkhorst-Pack grid that includes
         `\Gamma`
       * ``[(k11,k12,k13),(k21,k22,k23),...]``: Explicit list in units of the
         reciprocal lattice vectors
       * ``kpts=3.5``: `\vec k`-point density as in 3.5 `\vec k`-points per
         Å\ `^{-1}`.
   smearing: tuple
       The smearing of occupation numbers.  Must be a tuple:

       * ``('Fermi-Dirac', width)``
       * ``('Gaussian', width)``
       * ``('Methfessel-Paxton', width, n)``, where `n` is the order
         (`n=0` is the same as ``'Gaussian'``)

       Lower-case names are also allowed.  The ``width`` parameter is
       given in eV units.
   charge: float
      Charge of the system in units of `|e|` (``charge=1`` means one
      electron has been removed).  Default is ``charge=0``.
   nbands: int
      Number of bands.  Each band can be occupied by two electrons.

Not all of the above arguments make sense for all of ASE's
calculators.  As an example, Gromacs will not accept DFT related
keywords such as ``xc`` and ``smearing``.  In addition to the keywords
mentioned above, each calculator may have native keywords that are
specific to only that calculator.

Keyword arguments can also be set or changed at a later stage using
the :meth:`set` method:

.. method:: set(key1=value1, key2=value2, ...)


.. toctree::

   eam
   emt
   abinit
   amber
   castep
   cp2k
   crystal
   demon
   demonnano
   dftb
   dmol
   espresso
   exciting
   FHI-aims
   fleur
   gamess_us
   gaussian
   gromacs
   gulp
   socketio/socketio
   jacapo
   kim
   lammps
   lammpsrun
   mopac
   nwchem
   octopus
   onetep
   openmx
   orca
   psi4
   qchem
   siesta
   turbomole
   vasp
   qmmm
   checkpointing
   mixing
   loggingcalc
   dftd3
   others
   test
   ase_qmmm_manyqm
   ace

.. _calculator interface:

Calculator interface
====================

All calculators must have the following interface:

.. autoclass:: ase.calculators.interface.Calculator
   :members:


Electronic structure calculators
================================

These calculators have wave functions, electron densities, eigenvalues
and many other quantities.  Therefore, it makes sense to have a set of
standard methods for accessing those quantities:

.. autoclass:: ase.calculators.interface.DFTCalculator
   :members:
