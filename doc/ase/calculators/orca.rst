.. module:: ase.calculators.orca

======
ORCA
======

`ORCA <https://orcaforum.kofo.mpg.de/app.php/portal>`_ is a computational chemistry code 
that can do SCF, (TD)DFT, semi-empirical potentials, MP2, CASSCF, Coupled Cluster
calculations, and more. 


It is closed source, but free for academic users. Register on the forum to receive 
a download link for the binaries, as well as access ot the latest manual.


Many input examples are available at the 
`ORCA Input Library <https://sites.google.com/site/orcainputlibrary>`_.


.. highlight:: none

The :class:`ORCA` ASE-interface is very simple. Two keywords are defined::

  orcasimpleinput: str
      What you'd put after the "!" in an orca input file.

  orcablock: str
      What you'd put in the "% ... end"-blocks.


The ASE-calculator also works with the 
:mod:`~ase.calculators.qmmm.EIQMMM`-calculator 
for QM/MM simulations (see :mod:`~ase.calculators.qmmm` for 
more info). 

Setup and usage
===============

.. highlight:: bash

The default command that ASE will use to start ORCA is
``orca PREFIX.inp > PREFIX.out``. 

You can change this command by setting the
environment variable :envvar:`$ASE_ORCA_COMMAND`. (For example, add a line
to your ``.bashrc`` with ``export ASE_ORCA_COMMAND="my new command"``). 
This can be useful since the parallel MPI version of orca can require the full
path to the executable to be specified. 

.. highlight:: python

Orca wants to decide which sub-processes to parallelize via MPI itself, so you'd
almost always want a string in your ``orcablocks`` specifying the number of 
cores for the simulation, e.g.::

  from ase.calculators.orca import ORCA

  calc = ORCA(label='orcacalc', 
              orcasimpleinput='B3LYP def2-TZVP'
              orcablocks='%pal nprocs 16 end'

for a B3LYP/def2-TZVP calculation on 16 cores. 

