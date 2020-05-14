.. module:: ase.calculators.siesta

======
SIESTA
======

Introduction
============

SIESTA_ is a density-functional theory code for very large systems
based on atomic orbital (LCAO) basis sets.


.. _SIESTA: https://departments.icmab.es/leem/siesta/



Environment variables
=====================

The environment variable :envvar:`ASE_SIESTA_COMMAND` must hold the command
to invoke the siesta calculation. The variable must be a string 
where ``PREFIX.fdf``/``PREFIX.out`` are the placeholders for the 
input/output files. This variable allows you to specify serial or parallel 
execution of SIESTA.
Examples: ``siesta < PREFIX.fdf > PREFIX.out`` and
``mpirun -np 4 /bin/siesta4.0 < PREFIX.fdf > PREFIX.out``.

A default directory holding pseudopotential files :file:`.vps/.psf` can be
defined to avoid defining this every time the calculator is used.
This directory can be set by the environment variable
:envvar:`SIESTA_PP_PATH`.

Set both environment variables in your shell configuration file:

.. highlight:: bash

::

  $ export ASE_SIESTA_COMMAND="siesta < PREFIX.fdf > PREFIX.out"
  $ export SIESTA_PP_PATH=$HOME/mypps

.. highlight:: python

Alternatively, the path to the pseudopotentials can be given in
the calculator initialization. Listed below all parameters 
related to pseudopotential control.

===================== ========= ============= =====================================
keyword               type      default value description
===================== ========= ============= =====================================
``pseudo_path``       ``str``   ``None``      Directory for pseudopotentials to use
                                              None means using $SIESTA_PP_PATH
``pseudo_qualifier``  ``str``   ``None``      String for picking out specific type
                                              type of pseudopotentials. Giving
                                              ``example`` means that
                                              ``H.example.psf`` or
                                              ``H.example.vps`` will be used. None
                                              means that the XC.functional keyword
                                              is used, e.g. ``H.lda.psf``
``symlink_pseudos``   ``bool``  ``---``       Whether psedos will be sym-linked 
                                              into the execution directory. If 
                                              False they will be copied in stead.
                                              Default is True on Unix and False on
                                              Windows.
===================== ========= ============= =====================================


SIESTA Calculator
=================

These parameters are set explicitly and overrides the native values if different.

================ ========= =================== =====================================
keyword          type      default value       description
================ ========= =================== =====================================
``label``        ``str``   ``'siesta'``        Name of the output file
``mesh_cutoff``  ``float`` ``200*Ry``          Mesh cut-off energy in eV
``xc``           ``str``   ``'LDA'``           Exchange-correlation functional.
                                               Corresponds to either XC.functional
                                               or XC.authors keyword in SIESTA
``energy_shift`` ``float`` ``100 meV``         Energy shift for determining cutoff
                                               radii
``kpts``         ``list``  ``[1,1,1]``         Monkhorst-Pack k-point sampling
``basis_set``    ``str``   ``DZP``             Type of basis set ('SZ', 'DZ', 'SZP',
                                               'DZP')
``spin``         ``str``   ``'non-polarized'`` The spin approximation used, must be
                                               either ``'non-polarized'``, 
                                               ``'collinear'``, ``'non-collinear'``
                                               or ``'spin-orbit'``.
``species``      ``list``  ``[]``              A method for specifying the basis set  
                                               for some atoms.
================ ========= =================== =====================================

Most other parameters are set to the default values of the native interface.

Extra FDF parameters
====================

The SIESTA code reads the input parameters for any calculation from a
:file:`.fdf` file. This means that you can set parameters by manually setting
entries in this input :file:`.fdf` file. This is done by the argument:

>>> Siesta(fdf_arguments={'variable_name': value, 'other_name': other_value})

For example, the ``DM.MixingWeight`` can be set using

>>> Siesta(fdf_arguments={'DM.MixingWeight': 0.01})

The explicit fdf arguments will always override those given by other
keywords, even if it breaks calculator functionality.
The complete list of the FDF entries can be found in the official `SIESTA
manual`_.

.. _SIESTA manual: https://departments.icmab.es/leem/siesta/Documentation/Manuals/manuals.html

Example
=======

Here is an example of how to calculate the total energy for bulk Silicon,
using a double-zeta basis generated by specifying a given energy-shift:

>>> from ase import Atoms
>>> from ase.calculators.siesta import Siesta
>>> from ase.units import Ry
>>>
>>> a0 = 5.43
>>> bulk = Atoms('Si2', [(0, 0, 0),
...                      (0.25, 0.25, 0.25)],
...              pbc=True)
>>> b = a0 / 2
>>> bulk.set_cell([(0, b, b),
...                (b, 0, b),
...                (b, b, 0)], scale_atoms=True)
>>>
>>> calc = Siesta(label='Si',
...               xc='PBE',
...               mesh_cutoff=200 * Ry,
...               energy_shift=0.01 * Ry,
...               basis_set='DZ',
...               kpts=[10, 10, 10],
...               fdf_arguments={'DM.MixingWeight': 0.1,
...                              'MaxSCFIterations': 100},
...               )
>>> bulk.calc = calc
>>> e = bulk.get_potential_energy()

Here, the only input information on the basis set is, that it should
be double-zeta (``basis='DZP'``) and that the confinement potential
should result in an energy shift of 0.01 Rydberg (the
``energy_shift=0.01 * Ry`` keyword). Sometimes it can be necessary to specify
more information on the basis set.

Defining Custom Species
=======================
Standard basis sets can be set by the keyword ``basis_set`` directly, but for
anything more complex than one standard basis size for all elements,
a list of ``species`` must be defined. Each specie is identified by atomic
element and the tag set on the atom.

For instance if we wish to investigate a H2 molecule and put a ghost atom
(the basis set corresponding to an atom but without the actual atom) in the middle
with a special type of basis set you would write:

>>> from ase.calculators.siesta.parameters import Specie, PAOBasisBlock
>>> from ase import Atoms
>>> from ase.calculators.siesta import Siesta
>>> atoms = Atoms(
...     '3H',
...     [(0.0, 0.0, 0.0),
...      (0.0, 0.0, 0.5),
...      (0.0, 0.0, 1.0)],
...     cell=[10, 10, 10])
>>> atoms.set_tags([0, 1, 0])
>>>
>>> basis_set = PAOBasisBlock(
... """1
... 0  2 S 0.2
... 0.0 0.0""")
>>>
>>> siesta = Siesta(
...     species=[
...         Specie(symbol='H', tag=None, basis_set='SZ'),
...         Specie(symbol='H', tag=1, basis_set=basis_set, ghost=True)])
>>>
>>> atoms.calc = siesta

When more species are defined, species defined with a tag has the highest priority.
General species with ``tag=None`` has a lower priority.
Finally, if no species apply
to an atom, the general calculator keywords are used.


Pseudopotentials
================

Pseudopotential files in the ``.psf`` or ``.vps`` formats are needed.
Pseudopotentials generated from the ABINIT code and converted to
the SIESTA format are available in the `SIESTA`_ website.
A database of user contributed pseudopotentials is also available there.

Optimized GGA–PBE pseudos and DZP basis sets for some common elements
are also available from the `SIMUNE`_ website.

You can also find an on-line pseudopotential generator_ from the
OCTOPUS code.

.. _SIMUNE: https://www.simuneatomistics.com/siesta-pro/siesta-pseudos-and-basis-database/

.. _generator: http://www.tddft.org/programs/octopus/wiki/index.php/Pseudopotentials


Species can also be used to specify pseudopotentials:

>>> specie = Specie(symbol='H', tag=1, pseudopotential='H.example.psf')

When specifying the pseudopotential in this manner, both absolute
and relative paths can be given.
Relative paths are interpreted as relative to the set 
pseudopotential path.

Restarting from an old Calculation
==================================

If you want to rerun an old SIESTA calculation, whether made using the ASE
interface or not, you can set the keyword ``restart`` to the siesta ``.XV``
file. The keyword ``ignore_bad_restart`` (True/False) will decide whether
a broken file will result in an error(False) or the whether the calculator
will simply continue without the restart file.

Choosing the coordinate format
==============================
If you are mainly using ASE to generate SIESTA files for relaxation with native
SIESTA relaxation, you may want to write the coordinates in the Z-matrix format
which will best allow you to use the advanced constraints present in SIESTA.

======================= ========= ============= =====================================
keyword                 type      default value description
======================= ========= ============= =====================================
``atomic_coord_format`` ``str``   ``'xyz'``     Choose between ``'xyz'`` and 
                                                ``'zmatrix'`` for the format that 
                                                coordinates will be written in.
======================= ========= ============= =====================================

TDDFT Calculations
==================

It is possible to run Time Dependent Density Functional Theory (TDDFT) using the 
`PYSCF-NAO <https://github.com/cfm-mpc/pyscf/tree/nao/pyscf/lib/nao>`_ code together 
with the SIESTA code. This code allows to run TDDFT up to 
thousand atoms with small computational ressources. Visit the 
`github <https://github.com/cfm-mpc/pyscf/tree/nao>`_ webpage for 
further informations about PYSCF-NAO.

Example of code to calculate polarizability of Na8 cluster,::

  from ase.units import Ry, eV, Ha
  from ase.calculators.siesta import Siesta
  from ase import Atoms
  import numpy as np
  import matplotlib.pyplot as plt

  # Define the systems
  Na8 = Atoms('Na8',
               positions=[[-1.90503810, 1.56107288, 0.00000000],
                          [1.90503810, 1.56107288, 0.00000000],
                          [1.90503810, -1.56107288, 0.00000000],
                          [-1.90503810, -1.56107288, 0.00000000],
                          [0.00000000, 0.00000000, 2.08495836],
                          [0.00000000, 0.00000000, -2.08495836],
                          [0.00000000, 3.22798122, 2.08495836],
                          [0.00000000, 3.22798122, -2.08495836]],
               cell=[20, 20, 20])

  # Siesta input
  siesta = Siesta(
              mesh_cutoff=150 * Ry,
              basis_set='DZP',
              pseudo_qualifier='',
              energy_shift=(10 * 10**-3) * eV,
              fdf_arguments={
                  'SCFMustConverge': False,
                  'COOP.Write': True,
                  'WriteDenchar': True,
                  'PAO.BasisType': 'split',
                  'SCF.DM.Tolerance': 1e-4,
                  'DM.MixingWeight': 0.01,
                  'MaxSCFIterations': 300,
                  'DM.NumberPulay': 4,
                  'XML.Write': True})

  Na8.calc = siesta
  e = Na8.get_potential_energy()
  freq, pol = siesta.get_polarizability_pyscf_inter(label="siesta",
                                                    jcutoff=7,
                                                    iter_broadening=0.15/Ha,
                                                    xc_code='LDA,PZ',
                                                    tol_loc=1e-6,
                                                    tol_biloc=1e-7,
                                                    freq = np.arange(0.0, 5.0, 0.05))
  # plot polarizability
  plt.plot(freq, pol[:, 0, 0].imag)
  plt.show()

Remark: 
-------

The PYSCF-NAO code is still under active development and to have access to
it with ASE you will need to use this PYSCF `fork <https://github.com/cfm-mpc/pyscf>`_ 
and use the branch nao. To summarize::

  git clone https://github.com/cfm-mpc/pyscf
  git fetch
  git checkout nao

Then you can follow the instruction of the `README <https://github.com/cfm-mpc/pyscf/blob/nao/pyscf/lib/nao/README.md>`_.
The installation is relatively easy, go to the lib directory::
  
  cd pyscf/pyscf/lib
  cp cmake_arch_config/cmake.arch.inc-your-config cmake.arch.inc
  mkdir build
  cd build
  cmake ..
  make

Then you need to add the pyscf directory to your PYTHONPATH

.. code-block:: none

  export PYTHONPATH=/PATH-TO-PYSCF/pyscf:$PYTHONPATH

.. _Siesta Raman:

Raman Calculations with SIESTA and PYSCF-NAO
============================================

It is possible to calculate the Raman spectra with SIESTA, PYSCF-NAO anf the
vibration module from ASE. Example with CO2,::

  from ase.units import Ry, eV, Ha
  from ase.calculators.siesta import Siesta
  from ase.calculators.siesta.siesta_raman import SiestaRaman
  from ase import Atoms
  import numpy as np

  # Define the systems
  # example of Raman calculation for CO2 molecule,
  # comparison with QE calculation can be done from
  # https://github.com/maxhutch/quantum-espresso/blob/master/PHonon/examples/example15/README

  CO2 = Atoms('CO2',
              positions=[[-0.009026, -0.020241, 0.026760],
                         [1.167544, 0.012723, 0.071808],
                         [-1.185592, -0.053316, -0.017945]],
              cell=[20, 20, 20])

  # enter siesta input
  # To perform good vibrational calculations it is strongly advised
  # to relax correctly the molecule geometry before to actually run the
  # calculations. Then to use a large mesh_cutoff and to have the option
  # PAO.SoftDefault turned on
  siesta = Siesta(
      mesh_cutoff=450 * Ry,
      basis_set='DZP',
      xc="GGA",
      pseudo_qualifier='gga',
      energy_shift=(10 * 10**-3) * eV,
      fdf_arguments={
          'SCFMustConverge': False,
          'COOP.Write': True,
          'WriteDenchar': True,
          'PAO.BasisType': 'split',
          "PAO.SoftDefault": True,
          'SCF.DM.Tolerance': 1e-4,
          'DM.MixingWeight': 0.01,
          'MaxSCFIterations': 300,
          'DM.NumberPulay': 4,
          'XML.Write': True,
          'DM.UseSaveDM': True})

  CO2.calc = siesta

  ram = SiestaRaman(CO2, siesta, nfree=4, label="siesta", jcutoff=7, iter_broadening=0.15/Ha,
          xc_code='LDA,PZ', tol_loc=1e-6, tol_biloc=1e-7, freq = np.arange(0.0, 5.0, 0.05))

  ram.run()
  ram.summary(intensity_unit_ram='A^4 amu^-1')
  ram.write_spectra(start=200, intensity_unit_ram='A^4 amu^-1')


Further Examples
================
See also ``ase/test/calculators/siesta/test_scripts`` for further examples
on how the calculator can be used.


Siesta Calculator Class
=======================

.. autoclass:: ase.calculators.siesta.siesta.Siesta
