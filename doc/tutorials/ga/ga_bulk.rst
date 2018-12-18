.. _ga_bulk_tutorial:

=================================
 GA crystal structure prediction
=================================

Here we will use the GA to predict a crystal structure of a given chemical composition.

The implementation is based on:

   | M. Van den Bossche, H. Grönbeck, and B. Hammer
   | `Tight-Binding Approximation-Enhanced Global Optimization`__
   | J. Chem. Theory Comput. 2018, 14, 2797−2807

   __ http://doi.org/10.1021/acs.jctc.8b00039
   
and has much the same functionality as USPEX by Oganov and coworkers.


Initial population
==================

.. literalinclude:: ga_bulk_start.py

Run the GA search
=================

Now we run the search



All the operators
-----------------

.. autoclass:: ase.ga.bulk_mutations.StrainMutation
.. autoclass:: ase.ga.bulk_mutations.SoftMutation
              
                  
