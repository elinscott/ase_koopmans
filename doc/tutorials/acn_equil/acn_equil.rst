.. _Acetonitrile MD Box Equillibration:

Equilibrating an MD box of acetonitrile
=======================================

In this tutorial we see how to perform a thermal equilibration of an MD
box of classical acetonitrile molecules using the Langevin module and 
the implementation of an acetonitrile force field in ASE.  

The procedure we use for the equilibration closely follows the one 
presented in the tutorial :ref:`TIPnP Water Box Equillibration`.  

Since the TIPnP type water interpotentials are for rigid
molecules, there are no intramolecular force terms, and we need to
constrain all internal degrees of freedom. For this, we're
using the RATTLE-type constraints of the :ref:`FixBondLengths` class to
constrain all internal atomic distances (O-H1, O-H2, and H1-H2) for
each molecule.

.. literalinclude:: acn_equil.py

.. note::

  The temperature calculated by ASE is assuming all degrees of freedom
  are available to the system. Since the constraints have removed the 3
  vibrational modes from each water, the shown temperature will be 2/3
  of the actual value.

The procedure for the TIP4P force field is the same, with the following
exception: the atomic sequence **must** be OHH, OHH, ... .

So to perform the same task using TIP4P, you simply have to import
that calculator instead:

::

    from ase.calculators.tip4p import TIP4P, rOH, angleHOH

More info about the TIP4P potential: :mod:`ase.calculators.tip4p`
