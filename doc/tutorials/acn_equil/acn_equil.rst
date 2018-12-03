.. _Acetonitrile MD Box Equilibration:

Equilibrating an MD box of acetonitrile
=======================================

In this tutorial we see how to perform a thermal equilibration of an MD
box of classical acetonitrile molecules using the Langevin module and 
the implementation of an acetonitrile force field in ASE.  

The acetonitrile force field implemented in ASE (:mod:`ase.calculators.acn`)
is an interaction potential between three-site linear molecules, in which 
the atoms of the methyl group are treated as a single site centered on the 
methyl carbon, i.e. hydrogens are not considered explicitely. For this reason, 
while setting up a box of acetonitrile one has to assign the mass of a methyl 
to the outer carbon atom. Furthermore, the calculator requires the atomic 
sequence to be MeCN ... MeCN or NCMeNCMe ... NCMe, where Me represents the 
methyl site.

As for the TIPnP models, the acetonitrile potential works with rigid molecules.
However, due to the linearity of the acetonitrile molecular model, we cannot 
fix the geometry by constraining all interatomic distances using :class:`FixBondLengths`, 
as is done for TIPnP water. Instead, we must use the class :class:`FixLinearTriatomic` 


The MD procedure we use for the equilibration closely follows the one 
presented in the tutorial :ref:`TIPnP Water Box Equillibration`. A difference
is represented by the specification of the keyword *selectlinear*, which is
used to select molecules in the simulation box that are not linear triatomics. 
In this case, since all acetonitrile molecules are linear triatomics *selectlinear*
is given the value of an empty list.  

.. literalinclude:: acn_equil.py
