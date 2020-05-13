.. module:: ase.calculators.gaussian

========
Gaussian
========

`Gaussian <http://gaussian.com>`_ is a computational chemistry code
based on gaussian basis functions.


Setup
=====

.. highlight:: bash

Ask your system administrator to install Gaussian for you.

The ASE Gaussian calculator has been written with Gaussian 16 (g16) in mind,
but it will likely work with newer and older versions of Gaussian as well.
By default, the Calculator will look for executables named ``g16``, ``g09``,
and ``g03`` in that order. If your Gaussian executable is named differently,
or if it is not present in :envvar:`PATH`, then you must pass the path and name
of your Gaussian executable to the ``command`` keyword argument of the Gaussian
calculator. The default command looks like ``g16 < PREFIX.com > PREFIX.log``,
so template the ``command`` similarly. Alternatively, you may
set the :envvar:`ASE_GAUSSIAN_COMMAND` environment variable to the full
Gaussian executable command.


Examples
========

Here is a command line example of how to optimize the geometry of a
water molecule using the PBE density functional::

    $ ase build H2O | ase run gaussian -p xc=PBE,basis=3-21G -f 0.02
    $ ase gui stdin.traj@-1 -tg "a(1,0,2),d(0,1)"
    102.58928991353669 1.0079430292939233

.. highlight:: python

An example of creating a Gaussian calculator in the python interface is::

  from ase.calculators.gaussian import Gaussian

  calc = Gaussian(label='calc/gaussian',
                  xc='B3LYP',
                  basis='6-31+G*',
                  scf='maxcycle=100')

Parameters
==========

.. highlight:: none

The Gaussian calculator has three main types of parameters:

1. `Link0 keywords <https://gaussian.com/link0/>`_
2. `Route section keywords <https://gaussian.com/route/>`_
3. ASE-specific keywords, or convenience keywords.

The Gaussian calculator maintains a list of Link0 keywords and ASE-specific
keywords. Any keyword not on one of those two lists is assumed to be a route
section keyword, and will be placed in the Gaussian input file accordingly.

For example, consider the following Gaussian input file::

  %mem=1GB
  %chk=MyJob.chk
  %save
  #P b3lyp/6-31G scf=qc

  My job label

  0 1
  H 0.00 0.00 0.00
  H 0.00 0.00 0.74
  

.. highlight:: python

This would be generated with the following Python code::

  from ase import Atoms
  from ase.calculators.gaussian import Gaussian

  atoms = Atoms('H2', [[0, 0, 0], [0, 0, 0.74]])
  atoms.calc = Gaussian(mem='1GB',
                        chk='MyJob.chk',
                        save=None,
                        method='b3lyp',
                        basis='6-31G',
                        scf='qc')
  atoms.get_potential_energy()

Alternatively, you may use the ``xc`` keyword in place of the ``method``
keyword. ``xc`` is almost identical to ``method``, except that ``xc`` can
translate between the common definitions of some exchange-correlation
functionals and Gaussian's name for those functions, for example PBE to PBEPBE.
The ``method`` keyword will not do any translation, whatever value you provide
to ``method`` will be written to the input file verbatim. If both are provided,
``method`` overrides ``xc``.

Note that the Gaussian calculator puts each route keyword on its own line,
though this should not affect the result of the calculation.

When a route section keyword has multiple arguments, it is usually written
like ``scf(qc,maxcycle=1000)`` in the Gaussian input file. There are at least
two ways of generating this with the Gaussian calculator:
``Gaussian(scf="qc,maxcycle=100")`` and
``Gaussian(scf=['qc', 'maxcycle=100'])``, with the latter being somewhat more
convenient for scripting purposes.

Aside from the link-line and route section arguments, the Gaussian calculator
accepts a few additional convenience arguments.

============== ======== =============== ==================================================
keyword        type     default value   description
============== ======== =============== ==================================================
``label``      ``str``  ``'Gaussian'``  Name to use for input and output files.
``method``     ``str``  None            Level of theory to use, e.g. ``hf``, ``ccsd``,
                                        ``mp2``, or ``b3lyp``.  Overrides ``xc``
                                        (see below).
``xc``         ``str``  None            Level of theory to use. Translates several XC
                                        functionals from their common name (e.g. PBE) to
                                        their internal Gaussian name (e.g. PBEPBE).
``basis``      ``str``  None            The basis set to use. If not provided, no basis
                                        set will be requested, which usually results
                                        in STO-3G.  Maybe omitted if ``basisfile`` is set
                                        (see below).
``charge``     ``int``  See description The system charge. If not provided, it will be
                                        automatically determined from the Atoms object's
                                        ``initial_charges``.
``mult``       ``int``  See description The system multiplicity (spin + 1). If not
                                        provided, it will be automatically determined from
                                        the Atoms object's ``initial_magnetic_moments``.
``basisfile``  ``str``  None            The basis file to use. If a value is provided,
                                        ``basis`` may be omitted (it will be automatically
                                        set to ``'gen'``)
``extra``      ``str``  None            Extra lines to be included in the route section
                                        verbatim. It should not be necessary to use this,
                                        but it is included for backwards compatibility.
``addsec``     ``str``  None            Text to be added after the molecular geometry
                                        specification, e.g. for defining constraints
                                        with ``opt='modredundant'``.
``ioplist``    ``list`` None            A collection of IOPs definitions to be included in
                                        the route line.
============== ======== =============== ==================================================


GaussianOptimizer and GaussianIRC
=================================

There are also two Gaussian-specific :mod:`Optimizer <ase.optimize>`-like classes:
``GaussianOptimizer`` and ``GaussianIRC``, which can be used for geometry
optimizations and IRC calculations, respectively. These can be invoked in the
following way::

  from ase.calculators.gaussian import Gaussian, GaussianOptimizer
  atoms = ...
  calc_opt = Gaussian(...)
  opt = GaussianOptimizer(atoms, calc_opt)
  opt.run(fmax='tight', steps=100)

Note that this differs from ASE's standard Optimizer classes in a few key ways:

1. The ``fmax`` keyword takes a string rather than a force/energy criterion.
   Valid keywords are described in the
   `Gaussian manual page for optimization <http://gaussian.com/opt>`_.
2. Unlike ASE's standard Optimizer classes, it is not possible to iterate
   over the optimization with ``opt.irun(...)``.
3. It is also not possible to create a Trajectory file which records the
   optimization with ``opt = GaussianOptimizer(..., trajectory='opt.traj')``.
   However, it should be possible to obtain the trajectory by reading
   the Gaussian output file after the optimization has finished.

Additional arguments to Gaussian's ``opt`` keyword can be passed to the calculator
in the following way::

  opt.run(fmax='tight', steps=100, opt='calcfc,ts')

This example requests a Hessian calculation followed by optimization to a saddle point
("transition state optimization").

The ``GaussianIRC`` class can also be used to run IRC or pseudo-IRC calculations.
For example, the following script optimizes to a saddle point, then runs an IRC
optimization in the forward- and reverse-direction::

  from ase.calculators.gaussian import Gaussian, GaussianOptimizer, GaussianIRC
  atoms = ...

  # Optimize to a saddle point
  calc_opt = Gaussian(label='opt', ...)
  opt = GaussianOptimizer(atoms, calc_opt)
  opt.run(fmax='tight', steps=100, opt='calcfc,ts')
  tspos = atoms.positions.copy()

  # Do a vibrational frequency calculation and store the Hessian in a
  # checkpoint file, for use in subsequent IRC calculations
  atoms.calc = Gaussian(label='sp', chk='sp.chk', freq='')
  atoms.get_potential_energy()

  # Perform IRC in the "forwards" direction
  calc_irc_for = Gaussian(label='irc_for', chk='irc_for.chk', oldchk='sp.chk', ...)
  irc_for = GaussianIRC(atoms, calc_irc_for)
  irc_for.run(direction='forward', steps=20, irc='rcfc')  # reuses Hessian

  # Perform IRC in the "reverse" direction
  # First, restore TS positions
  atoms.positions[:] = tspos
  calc_irc_rev = Gaussian(label='irc_rev', chk='irc_rev.chk', oldchk='sp.chk', ...)
  irc_rev = GaussianIRC(atoms, calc_irc_rev)
  irc_rev.run(direction='reverse', steps=20, irc='rcfc')

It should also be possible to use the same ``Gaussian`` calculator object for
each of these steps, so long as the label is changed between calculations
(to avoid overwriting the output file) and the settings are changed appropriately.
It should also be possible to use the same ``GaussianIRC`` object for both the
forwards and reverse IRC calculations, so long as the label is changed (again to
avoid overwriting the output file).

.. autoclass:: Gaussian
