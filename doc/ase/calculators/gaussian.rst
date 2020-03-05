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
to your Gaussian executable to the ``command`` keyword argument of the Gaussian
calculator. The default command looks like ``g16 < PREFIX.com > PREFIX.log``,
so template the ``command`` to appear similar to that. Alternatively, you may
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
section keyword, and will be printed appropriately.

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

Alternatively, instead of the ``method`` keyword, you can use the ``xc``
keyword. ``xc`` is almost identical to ``method``, except that ``xc`` can
translate between the common definitions of some exchange-correlation
functionals and Gaussian's name for those functions, for example PBE to PBEPBE.
The ``method`` keyword will not do any translation, whatever value you provide
to ``method`` will be written to the input file verbatim.

Note that the Gaussian calculator puts each route keyword on its own line,
though this should not affect the result of the calculation.

When a route section keyword has multiple arguments, it is usually written
like ``scf(qc,maxcycle=1000)`` in the Gaussian input file. There are at least
two ways of generating this with the Gaussian calculator:
``Gaussian(scf="qc,maxcycle=100")`` and
``Gaussian(scf=['qc', 'maxcycle=100'])``, with the latter being somewhat more
convenient for scripting the generation of keyword arguments.

Aside from the link-line and route section arguments, the Gaussian calculator
accepts a few additional convenience arguments.

============== ======== =============== ==================================================
keyword        type     default value   description
============== ======== =============== ==================================================
``label``      ``str``  ``'Gaussian'``  Used to name input and output files.
``method``     ``str``  None            Used to specify the level of theory to use in the
                                        calculation, e.g. ``hf``, ``ccsd``, ``mp2``,
                                        or ``b3lyp``. ``method`` overrides ``xc`` (see
                                        below) if both are provided to the calculator.
``xc``         ``str``  None            Used to specify the level of theory, as an
                                        alternative  to ``method``. Accepts the common
                                        name of several exchange-correlation functionals,
                                        such as PBE and TPSS, and will automatically
                                        convert to Gaussian's internal name for these
                                        functionals (PBEPBE and TPSSTPSS, respectively).
``basis``      ``str``  None            The basis set to use for the calculation. If not
                                        provided, then no basis function will be added
                                        to the link line. This will usually result in
                                        ``sto-3g`` being used, unless the ``method``
                                        or ``xc`` specified implies a particular basis
                                        set. If you wish to use a basis set file, set
                                        ``basis='gen'``, and provide the path to your
                                        basis set file to the ``basisfile`` keyword
                                        (see below).
``charge``     ``int``  See description The charge of the system. If not provided, the
                                        charge will be set to the rounded sum of the
                                        initial atomic charges for your Atoms object
                                        (this will usually set the charge to 0).
``mult``       ``int``  See description The multiplicity of the system (spin + 1).
                                        If not provided, the multiplicity will be set
                                        to the rounded sum of the initial atomic
                                        magnetic moments of your Atoms object, plus 1
                                        (this will usually set the multiplicity to 1).
``basisfile``  ``str``  None            The full filesystem path to the basis file to use
                                        in the calculation, for when ``basis='gen'``.
                                        The basis file will be read in and inserted
                                        verbatim into the Gaussian input file, unless
                                        ``basisfile`` is formatted as in the example
                                        ``basisfile='@my-basisfile.gbs/N'``, in which
                                        case the value provided to ``basisfile`` will
                                        be inserted into the input file, instead of
                                        the contents of the referred file.
``extra``      ``str``  None            Any string passed to ``extra`` will be included
                                        verbatim in the route section of the input file.
                                        It should not be necessary to use this keyword,
                                        but it has been included for backwards
                                        compatibility with previous iterations of the
                                        Gaussian calculator.
``addsec``     ``str``  None            A string or collection of strings to be printed
                                        after the molecular geometry specification.
                                        For example, this can be used along with
                                        ``opt='modredundant'`` to fix particular internal
                                        coordinates during an optimization.
``ioplist``    ``list`` None            A collection of strings to be included in the
                                        route line, delimited by ``IOP(`` and ``)``.
                                        This can be used to change internal settings
                                        to Gaussian routines.
============== ======== =============== ==================================================

.. autoclass:: Gaussian
