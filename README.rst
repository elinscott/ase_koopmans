Atomic Simulation Environment for Koopmans
==========================================

ASE is a set of tools and Python modules for setting up, manipulating,
running, visualizing and analyzing atomistic simulations.

Webpage: http://wiki.fysik.dtu.dk/ase

ASE for Koopmans is a modified version of ASE that contains modifications
specific to Koopmans functionals. Ultimately we hope that these changes
will be merged with the official version of ASE.

ASE-Koopmans was forked from ASE v3.20.0b1.

Webpage: https://koopmans-functionals.org/

Requirements
------------

* Python_ 3.6 or later
* NumPy_ (base N-dimensional array package)
* SciPy_ (library for scientific computing)

Optional:


* For ASE's GUI: Matplotlib_ (2D Plotting)
* tkinter (for ase.gui)
* Flask (for ase.db web-interface)


Installation
------------

Add ``~/ase`` to your $PYTHONPATH environment variable and add
``~/ase/bin`` to $PATH (assuming ``~/ase`` is where your ASE folder is).


Testing
-------

Please run the tests::

    $ ase test  # takes 1 min.

and send us the output if there are failing tests.


Contact
-------

For ASE
~~~~~~~

* Mailing list: ase-users_
* IRC_: #ase on freenode.net

Please send us bug-reports, patches, code, ideas and questions.

For Koopmans
~~~~~~~~~~~~

* Mailing list: koopmans-users_

Example
-------

Geometry optimization of hydrogen molecule with NWChem:

>>> from ase import Atoms
>>> from ase.optimize import BFGS
>>> from ase.calculators.nwchem import NWChem
>>> from ase.io import write
>>> h2 = Atoms('H2',
               positions=[[0, 0, 0],
                          [0, 0, 0.7]])
>>> h2.calc = NWChem(xc='PBE')
>>> opt = BFGS(h2, trajectory='h2.traj')
>>> opt.run(fmax=0.02)
BFGS:   0  19:10:49    -31.435229     2.2691
BFGS:   1  19:10:50    -31.490773     0.3740
BFGS:   2  19:10:50    -31.492791     0.0630
BFGS:   3  19:10:51    -31.492848     0.0023
>>> write('H2.xyz', h2)
>>> h2.get_potential_energy()  # ASE's units are eV and Ang
-31.492847800329216

This example requires NWChem to be installed.

::

    $ ase gui h2.traj


.. _Python: http://www.python.org/
.. _NumPy: http://docs.scipy.org/doc/numpy/reference/
.. _SciPy: http://docs.scipy.org/doc/scipy/reference/
.. _Matplotlib: http://matplotlib.org/
.. _ase-users: https://listserv.fysik.dtu.dk/mailman/listinfo/ase-users
.. _koopmans-users: https://groups.google.com/g/koopmans-users
.. _IRC: http://webchat.freenode.net/?randomnick=0&channels=ase
