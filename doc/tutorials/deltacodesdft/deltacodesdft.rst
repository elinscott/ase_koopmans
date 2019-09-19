.. _dcdft tut:

Calculating Delta-values
========================

In this tutorial we compare the equation-of-state (EOS) calculated for 7 FCC
metals using values from :class:`~ase.calculators.emt.EMT`, WIEN2k and
experiment. Each EOS is described by three parameters:

* volume per atom
* bulk-modulus
* pressure derivative of bulk-modulus

Differences between two EOS'es can be measured by a single `\Delta` value
defined as:

.. math::

    \sqrt{\frac{\int_{V_a}^{V_b}
                (E_1(V) - E_2(V))^2 dV}
          {V_b - V_a}},

where `E_n(V)` is the energy per atom as a function of volume.
The `\Delta` value can be calculated using the
:func:`ase.utils.deltacodesdft.delta` function:

.. autofunction:: ase.utils.deltacodesdft.delta

.. seealso::

    * Collection of ground-state elemental crystals: :ref:`dcdft`
    * Equation-of-state module: :mod:`ase.eos`

We get the WIEN2k and experimental numbers from the :ref:`dcdft` ASE-collection
and we calculate the EMT EOS using this script:

.. literalinclude:: calculate.py

And fit to a Birch-Murnaghan EOS:

.. literalinclude:: fit.py

Result for Pt:

.. image:: Pt.png

Volumes in Ang^3:

.. csv-table::
    :file: volume.csv

Bulk moduli in GPa:

.. csv-table::
    :file: B.csv

Pressure derivative of bulk-moduli:

.. csv-table::
    :file: Bp.csv

Now, we can calculate `\Delta` between EMT and WIEN2k for Pt:

>>> from ase.utils.deltacodesdft import delta
>>> from ase.units import kJ
>>> delta(15.08, 278.67 * 1e-24 * kJ, 5.31,
...       15.64, 248.71 * 1e-24 * kJ, 5.46)
0.03205389052984122

Here are all the values (in meV/atom) calculated with the script below:

.. csv-table::
    :file: delta.csv

.. literalinclude:: tables.py
