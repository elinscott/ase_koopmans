.. _ga_bulk_tutorial:

=================================
 GA crystal structure prediction
=================================

Here we will use the GA to predict a crystal structure of a given chemical composition.

The implementation is based on:

   | M. Van den Bossche, H. Grönbeck, and B. Hammer
   | `Tight-Binding Approximation-Enhanced Global Optimization`__
   | J. Chem. Theory Comput. 2018, 14, 2797−2807

   __ https://doi.org/10.1021/acs.jctc.8b00039

and has much the same functionality as e.g. the USPEX program from the Oganov
group and the XtalOpt code by the Zurek group.

What sets crystal structure searches apart from other global optimization
problems, is that typically the cell vectors need to be treated as additional
degrees of freedom (in addition to the atomic coordinates). This results in
the following 'bulk GA'-specific functionality:

* Generating random structures now also involves randomized unit cell vectors.

* Cut-and-splice pairing needs to be done using scaled coordinates and
  includes a random combination of the two parent unit cells.

* A 'strain' mutation applies a random deformation of the parent unit cell,
  similar to how the 'rattle' mutation acts on the atomic coordinates.

* In a 'permustrain' mutation, atoms from different elements are swapped
  in addition to the 'strain' mutation.

* A 'soft' mutation displaces the atoms along a low-frequency vibrational
  mode obtained via e.g. a simple harmonic bond model.


As an example, we will search for the most energetically stable :mol:`Ag_2_4`
polymorphs using an EMT potential. Note that, in general, the number of atoms
per unit cell should be chosen carefully, as one will only find crystal structures
where the stoichiometry of the primitive cell is a divisor of the chosen
stoichiometry. As the stress tensor is currently not implemented in the EMT
calculator in ASE, the `ASAP <https://wiki.fysik.dtu.dk/asap/asap>`_ code will
be used instead.


Initial population
==================

The following script (:download:`ga_bulk_start.py`) creates a :file:`gadb.db`
database containing 20 randomly generated initial structures.

.. literalinclude:: ga_bulk_start.py


Run the GA search
=================

Now we can start the actual search (:download:`ga_bulk_run.py`),
which should only take a few minutes to complete.
The relaxation function, which performs the variable-cell
local optimization, is imported from :download:`ga_bulk_relax.py`.

.. literalinclude:: ga_bulk_run.py


All the bulk operators
----------------------

.. autoclass:: ase.ga.bulk_crossovers.CutAndSplicePairing
.. autoclass:: ase.ga.bulk_mutations.StrainMutation
.. autoclass:: ase.ga.bulk_mutations.PermuStrainMutation
.. autoclass:: ase.ga.bulk_mutations.SoftMutation
.. autoclass:: ase.ga.bulk_mutations.RotationalMutation
.. autoclass:: ase.ga.bulk_mutations.RattleRotationalMutation


All the helper functions
------------------------

.. autoclass:: ase.ga.bulk_utilities.CellBounds
.. autoclass:: ase.ga.bulk_startgenerator.StartGenerator
.. autoclass:: ase.ga.ofp_comparator.OFPComparator
