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

    @calcprop()
    def energy(self):
        return self['energy']

    @calcprop()
    def free_energy(self):
        return self['free_energy']

    @calcprop(shape=('natoms', 3))
    def forces(self):
        return self['forces']

    @calcprop(shape=(6,))
    def stress(self):
        stress = np.array(self['stress'])
        # Maybe offer to convert to Voigt form for lazy calculators
        return stress

    @calcprop(int)
    def nbands(self):
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
