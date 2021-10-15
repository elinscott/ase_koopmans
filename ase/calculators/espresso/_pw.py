"""Quantum ESPRESSO Calculator

export ASE_ESPRESSO_COMMAND="/path/to/pw.x -in PREFIX.pwi > PREFIX.pwo"

Run pw.x jobs.
"""


import numpy as np
import warnings
from ase.dft.kpoints import BandPath
from ase.calculators.calculator import PropertyNotPresent
from ._espresso import EspressoParent, EspressoWithBandstructure, error_template, warn_template


class Espresso(EspressoWithBandstructure, EspressoParent):
    """
    During initialisation, al options for pw.x are copied verbatim to
    the input file, and put into the correct section. Use ``input_data``
    for parameters that are already in a dict, all other ``kwargs`` are
    passed as parameters.

    Accepts all the options for pw.x as given in the QE docs, plus some
    additional options:

    input_data: dict
        A flat or nested dictionary with input parameters for pw.x
    pseudopotentials: dict
        A filename for each atomic species, e.g.
        ``{'O': 'O.pbe-rrkjus.UPF', 'H': 'H.pbe-rrkjus.UPF'}``.
        A dummy name will be used if none are given.
    kspacing: float
        Generate a grid of k-points with this as the minimum distance,
        in A^-1 between them in reciprocal space. If set to None, kpts
        will be used instead.
    kpts: (int, int, int), dict, or BandPath
        If kpts is a tuple (or list) of 3 integers, it is interpreted
        as the dimensions of a Monkhorst-Pack grid.
        If ``kpts`` is set to ``None``, only the Γ-point will be included
        and QE will use routines optimized for Γ-point-only calculations.
        Compared to Γ-point-only calculations without this optimization
        (i.e. with ``kpts=(1, 1, 1)``), the memory and CPU requirements
        are typically reduced by half.
        If kpts is a dict, it will either be interpreted as a path
        in the Brillouin zone (*) if it contains the 'path' keyword,
        otherwise it is converted to a Monkhorst-Pack grid (**).
        (*) see ase.dft.kpoints.bandpath
        (**) see ase.calculators.calculator.kpts2sizeandoffsets
    koffset: (int, int, int)
        Offset of kpoints in each direction. Must be 0 (no offset) or
        1 (half grid offset). Setting to True is equivalent to (1, 1, 1).


    .. note::
       Set ``tprnfor=True`` and ``tstress=True`` to calculate forces and
       stresses.

    .. note::
       Band structure plots can be made as follows:


       1. Perform a regular self-consistent calculation,
          saving the wave functions at the end, as well as
          getting the Fermi energy:

          >>> input_data = {<your input data>}
          >>> calc = Espresso(input_data=input_data, ...)
          >>> atoms.calc = calc
          >>> atoms.get_potential_energy()
          >>> fermi_level = calc.get_fermi_level()

       2. Perform a non-self-consistent 'band structure' run
          after updating your input_data and kpts keywords:

          >>> input_data['control'].update({'calculation':'bands',
          >>>                               'restart_mode':'restart',
          >>>                               'verbosity':'high'})
          >>> calc.set(kpts={<your Brillouin zone path>},
          >>>          input_data=input_data)
          >>> calc.calculate(atoms)

       3. Make the plot using the BandStructure functionality,
          after setting the Fermi level to that of the prior
          self-consistent calculation:

          >>> bs = calc.band_structure()
          >>> bs.reference = fermi_energy
          >>> bs.plot()
    """

    ext_in = '.pwi'
    ext_out = '.pwo'
    implemented_properties = ['energy', 'forces', 'stress', 'magmoms']
    command = 'pw.x -in PREFIX.pwi > PREFIX.pwo'

    def read_results(self):
        super().read_results()

        if isinstance(self.parameters.kpts, BandPath):
            # Add the bandstructure to the results. This is very un-ASE-y and should eventually be replaced
            self.band_structure(vbm_to_zero=False)

    def get_fermi_level(self):
        if self.calc is None:
            raise PropertyNotPresent(error_template % 'Fermi level')
        return self.calc.get_fermi_level()

    def get_ibz_k_points(self):
        if self.calc is None:
            raise PropertyNotPresent(error_template % 'IBZ k-points')
        ibzkpts = self.calc.get_ibz_k_points()
        if ibzkpts is None:
            warnings.warn(warn_template % 'IBZ k-points')
        return ibzkpts

    def get_k_point_weights(self):
        if self.calc is None:
            raise PropertyNotPresent(error_template % 'K-point weights')
        k_point_weights = self.calc.get_k_point_weights()
        if k_point_weights is None:
            warnings.warn(warn_template % 'K-point weights')
        return k_point_weights

    def get_eigenvalues(self, **kwargs):
        if self.calc is None:
            raise PropertyNotPresent(error_template % 'Eigenvalues')
        eigenvalues = self.calc.get_eigenvalues(**kwargs)
        if eigenvalues is None:
            warnings.warn(warn_template % 'Eigenvalues')
        return eigenvalues

    def get_number_of_spins(self):
        if self.calc is None:
            raise PropertyNotPresent(error_template % 'Number of spins')
        nspins = self.calc.get_number_of_spins()
        if nspins is None:
            warnings.warn(warn_template % 'Number of spins')
        return nspins

    def eigenvalues_from_results(self):
        assert hasattr(self.calc, 'kpts'), 'Please call {0}.calculate() prior to calling {0}.band_structure'.format(
            self.__class__.__name__)

        return np.array([[k.eps_n for k in self.calc.kpts if k.s == i_spin]
                         for i_spin in range(self.calc.get_number_of_spins())])
