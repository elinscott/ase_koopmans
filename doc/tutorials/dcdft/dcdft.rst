.. _dcdft tut:

Calculating Delta-values
========================

The `\Delta` value is defined as:

.. math::

    \sqrt{\int_{V_a}^{V_b}
          (E_1(V) - E_2(V))^2 / (V_b - V_a)}

and can be calculated using the :func:`ase.utils.dcdft.delta` function:

.. autofunction:: ase.utils.dcdft.delta

.. seealso::

    * :ref:`dcdft`
    * :mod:`ase.eos`

Let us calculate some `\Delta` values for some metals using :class:`~ase.calculators.emt.EMT`:

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

>>> from ase.utils.dcdft import delta
>>> from ase.units import kJ
>>> delta(15.08, 278.67 * 1e-24 * kJ, 5.31,
...       15.64, 248.71 * 1e-24 * kJ, 5.46)
0.03205389052984122

Here are all the values calculated with the script below:

.. csv-table::
    :file: delta.csv

.. literalinclude:: tables.py
