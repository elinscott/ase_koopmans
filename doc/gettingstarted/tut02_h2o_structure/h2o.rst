Structure optimization: :mol:`H_2O`
===================================

Let's calculate the structure of the :mol:`H_2O` molecule.

.. admonition:: Exercise

   Create an :class:`~ase.Atoms` object representing an :mol:`H_2O`
   molecule by providing chemical symbols and a guess for the positions.
   Visualize it, making sure the molecule is V shaped.


.. admonition:: Exercise

   Run a self-consistent calculation of the approximate H2O molecule
   using GPAW.


Optimizers
----------

We will next want to optimize the geometry.
ASE provides :mod:`several optimization algorithms <ase.optimize>`
that can run on top of :class:`~ase.Atoms` equipped with a calculator::

  from ase.optimize import BFGS
  opt = BFGS(atoms, trajectory='opt.traj', logfile='opt.log')
  opt.run(fmax=0.05)

.. admonition:: Exercise

   Run a structure optimization, thus calculating the equilibrium
   geometry of :mol:`H_2O`.

The ``trajectory`` keyword above ensures that the trajectory of intermediate
geometries is written to :file:`opt.traj`.

.. admonition:: Exercise

   Visualize the output trajectory and play it as an animation.
   Use the mouse to drag a box around and select the three atoms —
   this will display the angles between them.
   What is H–O–H angle of :mol:`H_2O`?

As always in ASE, we can do things programmatically, too,
if we know the right incantations::

  from ase.io import read
  atoms = read('opt.traj')
  print(atoms.get_angle(0, 1, 2))
  print(atoms.get_angle(2, 0, 1))
  print(atoms.get_angle(1, 2, 0))

The documentation on the :class:`~ase.Atoms` object provides
a long list of methods.


G2 molecule dataset
-------------------

ASE knows many common molecules, so we did not really need to type in
all the molecular coordinates ourselves.  As luck would have it, the
:func:`ase.build.molecule` function does exactly what we need::

  from ase.build import molecule
  atoms = molecule('H2O', vacuum=3.0)

This function returns a molecule from the G2 test set, which is nice
if we remember the exact name of that molecule, in this case `'H2O'`.
In case we don't have all the molecule names memorized, we can work
with the G2 test set using the more general :mod:`ase.collections.g2`
module::

  from ase.collections import g2

  print(g2.names)  # These are the molecule names
  atoms = g2['CH3CH2OH']
  view(atoms)
  view(g2)  # View all 162 systems


Use another calculator
----------------------

We could equally well substitute
another calculator, often accessed through imports like ``from
ase.calculators.emt import EMT`` or ``from ase.calculators.aims import
Aims``.  For a list, see :mod:`ase.calculators` or run::

  $ ase info --calculators

For these tutorials we also have FHI-Aims installed.  Let's run the
same relaxation with FHI-Aims then.  But in the list above, Aims
(probably) wasn't listed as available.
We first need to tell ASE how to run Aims
-- a typical configuration step for many ASE calculators.
This means specifying 1)
the command used to run Aims, and 2) where to find information about
chemical species.  We can do this by setting environment variables in the
shell:

::

   $ export ASE_AIMS_COMMAND=aims.x
   $ export AIMS_SPECIES_DIR=/home/alumne/software/FHIaims/species_defaults/light

Now ``ase info --calculators`` should tell us that it thinks
Aims is installed as ``aims.x``.

However, if we open a new shell it will forget this.  And we don't want to
modify ``.bashrc`` on these computers.  Let's instead set these variables
in our Python script:

::

   import os
   os.environ['ASE_AIMS_COMMAND'] = 'aims.x'
   os.environ['AIMS_SPECIES_DIR'] = '/home/alumne/software/FHIaims/species_defaults/light'


.. admonition:: Exercise

  Run a structure optimization of :mol:`H_2O`
  using the FHI-:class:`~ase.calculators.aims.Aims` calculator.

To enable the calculation of forces, you will need ``compute_forces=True``.
Aims will want an explicitly given XC functional, so we put ``xc='LDA'``.
The ``xc`` keyword is supported by several ASE calculators to make it easier
to specify common XC functionals.

After running the calculation, some new files will be present.
ASE has generated :file:`control.in` and :file:`geometry.in`, then
ran FHI-aims on them, producing :file:`aims.out`.
Be sure to briefly inspect the files.
Being perfectionist and/or paranoid,
we of course want to be sure that the ASE interface
set the parameters the way we wanted them.

Most ASE calculators can be made to generate a file
without triggering a calculation using ``calc.write_input_file(atoms)``.
This is useful, say, if you want to generate the files now but run them
later, with or without ASE.

ASE knows many file formats.  :func:`ase.io.read` can read both the
input file and the output file, returning :class:`~ase.Atoms`.
These files can also be opened directly with the ASE GUI.

Note that by default, subsequent calculations will overwrite each other.
Hence the Aims input and output files correspond to the final step of the
structure relaxation.  The documentation on :mod:`ase.optimize`
will tell us that we can override this behaviour by adding an observer,
or using the even more flexible :meth:`ase.optimize.Dynamics.irun` method
to force different steps into different directories.

Appendix: Communication between calculators and codes
-----------------------------------------------------

What follows is not necessary knowledge for normal usage of ASE. Unless
you are interested in how to optimize the communication between ASE and
external calculators you may skip ahead.

Different calculators communicate with computational codes in different ways.
GPAW is written in Python, so ASE and GPAW run within the same process.
However FHI-aims is a separate program.  What the Aims calculator
does for us is to generate an input file, run FHI-aims, read the output,
and return the results.

We just ran a relaxation which involved multiple geometry steps.  Each
step, a new Aims process is started and later stopped.  This is
inefficient because the ground-state density and wavefunctions of one
step would be an excellent initial guess for the next step, lowering
the number of steps necessary to converge.
But these quantities are lost when the program terminates.  To get
the best performance in structure optimisations and dynamics, we need to
avoid this loss of efficiency.

Many ASE calculators support more advanced ways of communicating.
These calculators can communicate with persistent external processes
over pipes (:class:`~ase.calculators.lammpsrun.Lammpsrun`, :class:`~ase.calculators.cp2k.CP2K`) or sockets (:class:`~ase.calculators.siesta.Siesta`,
:class:`~ase.calculators.aims.Aims`, :class:`~ase.calculators.espresso.Espresso`),
or they can work within the same process
through direct library calls
(:class:`~ase.calculators.lammpslib.Lammpslib`, GPAW).


ASE can communicate with FHI-aims over sockets using the i-PI protocol (http://ipi-code.org/).  This is done by wrapping the calculator in a
:class:`ase.calculators.socketio.SocketIOCalculator`.  The socket
calculator will use the calculator it wraps to launch a calculation,
then run it.

The documentation on the socket I/O calculator already provides full examples,
so we only need minor adjustments to run them on our local machine.

.. admonition:: Optional exercise

   Based on our previous relaxation with FHI-aims, write a script
   which runs the same calculation using the
   :class:`ase.calculators.socketio.SocketIOCalculator`.

You can run :command:`time python3 myscript.py` to see how long time
the calculation takes in total.  How much of a speedup do you get from
running the relaxation over a socket?  INET sockets often have high
latency.  If you don't see much of a speedup, this is probably why.
In that case, try switching to a UNIX socket.

The socket I/O calculator automatically generated an input file and also
immediately launched the calculation.  Since it only launches the process
once, subsequent steps don't overwrite each other and we can find all the
intermediate steps in :file:`aims.out`.

Solutions
---------

GPAW optimisation:

.. literalinclude:: solution/optimise.py

FHI-aims optimisation:

.. literalinclude:: solution/optimise_aims.py


FHI-aims/socket-io optimisation:

.. literalinclude:: solution/optimise_aims_socketio.py
