.. module:: ase.calculators.kim

===
KIM
===

.. note::
   This package requires the *KIM API* package, which is `hosted on GitHub
   <https://github.com/openkim/kim-api>`__ and available through many binary package
   managers.  See `openkim.org/kim-api <https://openkim.org/kim-api>`_ for installation
   options.

.. note::
   This package requires the *kimpy* python package, which is `hosted on GitHub
   <https://github.com/openkim/kimpy>`__ and also made available through `PyPI
   <https://pypi.org/project/kimpy/>`_.

The `Open Knowledgebase of Interatomic Models (OpenKIM) <https://openkim.org>`_ is an
NSF-funded project aimed at providing easy access to standardized implementations of
classical interatomic potentials that can be used with a variety of molecular simulation
codes.  This ASE package allows one to easily use any potential archived in OpenKIM through
ASE.  In order to explain its structure, we must first describe the two different types
of interatomic models in KIM: [#kimmodels]_ [#typesofkimcontent]_

* **Portable Models (PMs)**
  These models can be used with any KIM API-compliant simulator, either directly or
  through their corresponding ASE calculator.  Portable models published on openkim.org
  are identified by the string '__MO_' in their name.

* **Simulator Models (SMs)**
  These models are essentially just wrappers around a set of commands (and often one or
  more parameter files) in a specific simulator that are used to define the model.
  Simulator models published on openkim.org are identified by the string '__SM_' in their
  name.

These two types of KIM models require different calculators to work: PMs work through a
designated calculator that uses the `kimpy <https://github.com/openkim/kimpy>`_ library
(which provides a set of python bindings to the `KIM API <https://openkim.org/kim-api>`_)
in order to set up communication between ASE and the model.  This allows, for example,
the positions of the atoms and the neighbor list in ASE to be communicated to the model,
and for the energy and forces predicted by the model for that configuration to be
communicated back to ASE in a standard format.  On the other hand, SMs are
programmatically nothing more than a recording of some set of simulator commands and
usually one or more parameter files. [#smobject]_  Accordingly, they require supplying
these commands to the calculator(s) corresponding to the specific simulator associated
with the SM.

Because of this separation, the :mod:`ase.calculators.kim` package contains two modules:
:git:`ase.calculators.kim.kim <ase/calculators/kim/kim.py>` and
:git:`ase.calculators.kim.kimmodel <ase/calculators/kim/kimmodel.py>`. The first of these
contains a *wrapper function* named :func:`ase.calculators.kim.kim.KIM` that takes as
input the name of a KIM model, automatically determines whether it is a PM or an SM, then
constructs and returns an appropriate ASE calculator. [#getmodelsupportedspecies]_ For
example, if the name of an installed PM is passed, the ``KIM`` function will (by default)
initialize an instance of :class:`ase.calculators.kim.kimmodel.KIMModelCalculator` for it
and return it as its output.  If the name of a LAMMPS_-based SM is passed, the calculator
will (by default) return an instance of the :class:`ase.calculators.lammpslib.LAMMPSlib`
calculator.  The specific calculator type returned can be controlled using the
``simulator`` argument.

.. autofunction:: ase.calculators.kim.kim.KIM

========
Examples
========
By default, the KIM API will install several example models on your system, including a
portable model (PM) for argon named "ex_model_Ar_P_Morse_07C".  Suppose we wanted to know
the potential energy predicted by this model for an FCC argon lattice at a lattice
spacing of `a = 5.25`.  This can be accomplished like so:

::

    from ase.lattice.cubic import FaceCenteredCubic
    from ase.calculators.kim.kim import KIM

    atoms = FaceCenteredCubic(symbol='Ar', latticeconstant=5.25, size=(1,1,1))
    calc = KIM("ex_model_Ar_P_Morse_07C")
    atoms.set_calculator(calc)

    energy = atoms.get_potential_energy()
    print("Potential energy: {} eV".format(energy))

Because the ``simulator`` keyword argument is not specified, the object ``calc`` returned
here is actually an instance of :class:`ase.calculators.kim.kimmodel.KIMModelCalculator`
(because our model is a portable model) and uses the neighbor list library implemented in
kimpy.  If we wanted to use ASE's internal neighbor list mechanism, we could specify it
by modifying the corresponding line to:

::

    calc = KIM("ex_model_Ar_P_Morse_07C", options={"ase_neigh": True})

If, for some reason, we want to make sure that the ASE LAMMPS calculator
(:class:`ase.calculators.lammpsrun.LAMMPS`) is the specific one used to run our portable
model, we can specify it using the ``simulator`` argument:

::

    calc = KIM("ex_model_Ar_P_Morse_07C", simulator="lammpsrun")

Using a KIM simulator model requires no additional effort.  Using the example
LAMMPS-based simulator model bundled with the KIM API,
"Sim_LAMMPS_LJcut_AkersonElliott_Alchemy_PbAu":

::

    from ase.lattice.cubic import FaceCenteredCubic
    from ase.calculators.kim.kim import KIM

    atoms = FaceCenteredCubic(symbol='Au', latticeconstant=4.07, size=(1,1,1))
    calc = KIM("Sim_LAMMPS_LJcut_AkersonElliott_Alchemy_PbAu")
    atoms.set_calculator(calc)

    energy = atoms.get_potential_energy()
    print("Potential energy: {} eV".format(energy))

In this case, because ``simulator`` was not specified, the default behavior is that the
object ``calc`` returned is an instance of :class:`ase.calculators.lammpslib.LAMMPSlib`.

.. _LAMMPS: http://lammps.sandia.gov

.. rubric:: Footnotes

.. [#kimmodels] A listing of all models currently available in KIM can be found `here <https://openkim.org/browse/models/alphabetical>`__.

.. [#typesofkimcontent] See `here <https://openkim.org/doc/repository/kim-content/>`_ for more details about different types of content in KIM.

.. [#smobject] Simulator Models (SMs) are actually always used in the form of a binary shared object.  This shared object contains not only the commands and other metadata related to the SM within it, but also has any relevant parameter files (if any) directly compiled into it.

.. [#getmodelsupportedspecies] As a convenience, this module also contains the function :func:`ase.calculators.kim.kim.get_model_supported_species`, which retrieves a tuple of all of the species that a model can compute interactions for.
