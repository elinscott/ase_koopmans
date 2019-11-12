import functools

import numpy as np


all_properties = set()

def calcprop(dtype=float, shape=tuple()):
    def checker(func):
        all_properties.add(func.__name__)
        @functools.wraps(func)
        def prop(self):
            obj = func(self)
            if not shape:
                return dtype(obj)

            actual_shape = tuple(self[key] if isinstance(key, str) else key
                                 for key in shape)
            assert np.shape(obj) == actual_shape
            obj = np.array(obj)
            assert obj.shape == actual_shape
            assert obj.dtype == dtype
            return obj

            #assert np.issubddtype(dtype)
            return obj
        return property(prop)
    return checker


class CalculatorOutputs:
    def __init__(self, results):
        self.results = results

    def __getitem__(self, key):
        return self.results[key]

    def __contains__(self, key):
        return key in self.results

    def __iter__(self):
        return iter(self.results)

    @calcprop()
    def energy(self):
        return self['energy']

    @calcprop()
    def free_energy(self):
        return self['free_energy']

    @calcprop(shape=('natoms',))
    def energies(self):
        return self['energies']

    # energies?
    # (free_energies?  But that's unheard of.)

    @calcprop(shape=('natoms', 3))
    def forces(self):
        # We don't necessarily want to store the whole Atoms object
        # (after all our purpose is to represent the outputs!)
        return self['forces']

    @calcprop(shape=(6,))
    def stress(self):
        stress = np.array(self['stress'])
        # Maybe offer to convert to Voigt form for lazy calculators
        return stress

    @calcprop(int)
    def nbands(self):
        # We can be more intelligent here -- if we have self['eigenvalues'],
        # we necessarily have nbands.
        #
        # But we can also rely on a normalization step to guarantee
        # the presence of these things.
        return self['nbands']

    @calcprop(int)
    def nspins(self):
        return self['nspins']

    @calcprop(int)
    def nkpts(self):
        return self['nkpts']

    @calcprop(shape=('nspins', 'nkpts', 'nbands'))
    def eigenvalues(self):
        return self['eigenvalues']

    @calcprop(shape=('nspins', 'nkpts', 'nbands'))
    def occupations(self):
        return self['occupations']

    @calcprop()
    def fermi_level(self):
        return self['fermi_level']

    @calcprop(shape=('nkpts', 3))
    def ibz_kpoints(self):
        return self['ibz_kpoints']

    @calcprop(shape=('nkpts',))
    def kpoint_weights(self):
        return self['kpoint_weights']

    #def to_singlepoint(self, atoms):
    #    from ase.calculators.singlepoint import SinglePointDFTCalculator
    #    return SinglePointDFTCalculator(atoms,
    #                                    efermi=self.fermi_level,

    # We can also retrieve (P)DOS and band structure.  However:
    #
    # * Band structure may require bandpath, which is an input, and
    #   may not necessarily be easy or possible to reconstruct from
    #   the outputs.
    #
    # * Some calculators may produce the whole BandStructure object in
    #   one go (e.g. while parsing)
    #
    # * What about HOMO/LUMO?  Can be obtained from
    #   eigenvalues/occupations, but some codes provide real data.  We
    #   probably need to distinguish between HOMO/LUMO inferred by us
    #   versus values provided within the output.
    #
    # * HOMO is sometimes used as alternative reference energy for
    #   band structure.
    #
    # * What about spin-dependent (double) Fermi level?
    #
    # * What about 3D arrays?  We will almost certainly want to be
    #   connected to an object that can load dynamically from a file.
