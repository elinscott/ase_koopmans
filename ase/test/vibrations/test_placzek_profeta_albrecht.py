"""
Test resonant Raman implementations
"""

import numpy as np
from ase import Atoms
from ase.calculators.lj import LennardJones
from ase.vibrations.placzek import Placzek, Profeta
from ase.vibrations.albrecht import Albrecht
from ase.calculators.excitation import ExcitationList, Excitation


class FakeExcitation(Excitation):
    def __init__(self, energy, mur, muv, magn=None):
        self.energy = energy
        self.fij = 1.
        self.mur = mur
        self.muv = muv
        self.magn = magn

    def outstring(self):
        string = '{0:g}   '.format(self.energy)

        def format_me(me):
            string = ''
            if me.dtype == float:
                for m in me:
                    string += ' {0:.5e}'.format(m)
            else:
                for m in me:
                    string += ' {0.real:.5e}{0.imag:+.5e}j'.format(m)
            return string

        string += '  ' + format_me(self.mur)
        if self.muv is not None:
            string += '  ' + format_me(self.muv)
        if self.magn is not None:
            string += '  ' + format_me(self.magn)
        string += '\n'

        return string


class FakeExcited(ExcitationList):
    def __init__(self, calculator):
        self.calculator = calculator
        energies = np.arange(1, 4)
        murs = np.eye(3) / energies
        muvs = np.eye(3)
        for energy, mur, muv in zip(energies, murs, muvs):
            self.append(FakeExcitation(energy, mur, muv))
    def write(self, fname):
        with open(fname, 'w') as f:
            for ex in self:
                 f.write(ex.outstring())

lj = LennardJones(sigma=2**(1/6))
atoms = Atoms('H2', positions=[[0, 0 ,0], [0, 0, 1]])
atoms.set_calculator(lj)

pz = Placzek(atoms, FakeExcited)
pz.run()

om = 0.1
pzi = pz.absolute_intensity(omega=om)[-1]
print(pzi)
