Atoms and calculators
=====================

ASE allows atomistic calculations to be scripted with different
computational codes. In this introductory exercise, we go through the
basic concepts and workflow of ASE and will eventually
calculate the binding curve of :mol:`N_2`.

These tutorials often use the electronic structure code GPAW.
They can be completed just as well
using other supported codes, subject to minor adjustments.


Python
------

In ASE, calculations are performed by writing and running Python
scripts.  A very short primer on Python can be found in the
:ref:`ASE documentation <what is python>`.
If you are new to Python it would be wise to look through
this to understand the basic syntax, datatypes, and
things like imports.  Or you can just wing it --- we won't judge.


Atoms
-----

Let's set up a molecule and run a DFT calculation.
We can create
simple molecules by manually typing the chemical symbols and a
guess for the atomic positions in Ångström.  For example
:mol:`N_2`::

  from ase import Atoms
  atoms = Atoms('N2', positions=[[0, 0, -1], [0, 0, 1]])

Just in case we made a mistake, we should visualize our molecule
using the :mod:`ASE GUI <ase.gui>`::

  from ase.visualize import view
  view(atoms)

Equivalently we can save the atoms in some format, often ASE's own
:mod:`~ase.io.trajectory` format::

  from ase.io import write
  write('myatoms.traj', atoms)

Then run the GUI from a terminal::

  $ ase gui myatoms.traj

ASE supports quite a few different formats.   For the full list, run::

  $ ase info --formats

Although we won't be using all the ASE commands any time soon,
feel free to get an overview::

  $ ase --help

.. admonition:: Exercise

   Write a script which sets up and saves an :mol:`N_2` molecule,
   then visualize it.


Calculators
-----------

Next let us perform an electronic structure calculation.  ASE uses
:mod:`~ase.calculators` to perform calculations.  Calculators are
abstract interfaces to different backends which do the actual computation.
Normally, calculators work by calling an external electronic structure
code or force field code.  To run a calculation, we must first create a
calculator and then attach it to the :class:`~ase.Atoms` object. Here we
use GPAW and set a few calculation parameters as well::

  from gpaw import GPAW

  calc = GPAW(mode='lcao', basis='dzp', txt='gpaw.txt', xc='LDA')
  atoms.calc = calc


Different electronic structure codes have different input parameters.
GPAW can use real-space grids
(``mode='fd'``), planewaves (``mode='pw'``), or localized atomic orbitals
(``mode='lcao'``) to represent the wavefunctions.
Here we have asked for the faster but
less accurate LCAO mode, together with the standard double-zeta polarized basis
set (``'dzp'``).  GPAW and many other codes require a unit cell (or simulation
box) as well.  Hence we center the atoms within a box, leaving 3 Å
of empty space around each atom::

  atoms.center(vacuum=3.0)
  print(atoms)

The printout will show the simulation box (or ``cell``) coordinates,
and the box can also be seen in the GUI.

Once the :class:`~ase.Atoms` have a calculator with appropriate parameters,
we can do things like calculating energies and forces::

  e = atoms.get_potential_energy()
  print('Energy', e)
  f = atoms.get_forces()
  print('Forces')
  print(f)

This will give us the energy in eV and the forces in eV/Å.
(ASE also provides ``atoms.get_kinetic_energy()``, referring to the kinetic
energy of the nuclei if they are moving.  In DFT calculations,
we normally just want the Kohn--Sham ground state energy which is the
"potential" energy as provided by the calculator.)

Calling ``get_potential_energy()`` or ``get_forces()`` triggers a
selfconsistent calculation and gives us a lot of output text.
Inspect the :file:`gpaw.txt` file.  You can review the text file to see what
computational parameters were chosen.  Note how the ``get_forces()``
call did not actually trigger a *new* calculation --- the forces
were evaluated from the ground state that was already calculated,
so we only ran one calculation.


Binding curve
-------------

The strong point of ASE is that things are scriptable.
``atoms.positions`` is a numpy array with the atomic positions::

  print(atoms.positions)

We can move the atoms by adding or assigning other values into some of the
array elements.  Then we can trigger a new calculation by calling
``atoms.get_potential_energy()`` or
``atoms.get_forces()`` again.

.. admonition:: Exercise

   Move one of the atoms so you trigger two calculations in one script.

This way we can implement any series of calculations.  When running
multiple calculations, we often want to write them into a file.
We can use the standard trajectory format to write multiple calculations
(atoms and energy) like this::

  from ase.io.trajectory import Trajectory
  traj = Trajectory('mytrajectory.traj', 'w')
  ...
  traj.write(atoms)

.. admonition:: Exercise

   Write a loop, displacing one of the atoms in small steps to
   trace out a binding energy curve :math:`E(d)` around the equilibrium
   distance.  Save each step to a trajectory and visualize.
   What is the equilibrium distance?

Be careful that the atoms don't move too close to the edge of the
simulation box (or the electrons will squeeze against the box, increasing
energy and/or crashing the calculation).

.. note::

   The binding will be much too strong because our calculations are
   spin-paired, and the atoms would polarise as they move apart.

In case we forgot to write the trajectory,
we can also run ASE GUI on the :file:`gpaw.txt` file although its
printed precision is limited.

Although the GUI will plot the energy curve for us, publication
quality plots usually require some manual tinkering.
ASE provides two functions to read trajectories or other files:

 * :func:`ase.io.read` reads and returns the last image, or possibly a list of images if the ``index`` keyword is also specified.

 * :func:`ase.io.iread` reads multiple images, one at a time.

Use :func:`ase.io.iread` to read the images back in, e.g.::

  for atoms in iread('mytrajectory.traj'):
      print(atoms)

.. admonition:: Exercise

   Plot the binding curve (energy as a function of distance) with matplotlib.
   You will need to collect the energies and the distances when looping
   over the trajectory.  The atoms already have the energy.  Hence, calling
   ``atoms.get_potential_energy()`` will simply retrieve the energy
   without calculating anything.



.. admonition:: Optional exercise

   To get a more correct binding energy, set up an isolated N atom
   and calculate its energy.  Then calculate the molecular
   atomisation energy
   :math:`E_{\mathrm{atomisation}} = E[\mathrm N_2] - 2 E[\mathrm N]`
   of the :mol:`N_2` molecule.

   You can use ``atoms.set_initial_magnetic_moments([3])`` before
   triggering the calculation to tell
   GPAW that your atom is spin polarized.


Solutions
---------

.. literalinclude:: solution/binding_curve.py

.. literalinclude:: solution/plot_binding_curve.py

.. literalinclude:: solution/spinpol.py
