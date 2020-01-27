"""
Test resonant Raman implementations
"""
import sys
import numpy as np
from ase import Atoms
from ase.calculators.lj import LennardJones
from ase.vibrations import Vibrations

from ase.vibrations.placzek import Placzek, Profeta
from ase.vibrations.albrecht import Albrecht
from ase.calculators.excitation import ExcitationList, Excitation
from ase.calculators.h2lj import H2LJ, H2LJExcitedStates


def test_placzek_run():
    atoms = H2LJ()
    pz = Placzek(atoms, H2LJExcitedStates, txt='-')
    pz.run()


test_placzek_run()

sys.exit()


print('??? 2')

om = 0.1
pzi = pz.absolute_intensity(omega=om)[-1]
print(pzi)

sys.exit()

class FakeExcitation(Excitation):
    def __init__(self, energy=None, mur=None, muv=None, magn=None, string=None):
        if string is not None:
            self.fromstring(string)
        else:
            self.energy = energy
            self.mur = mur
            self.muv = muv
            self.magn = magn
        self.fij = 1.
        self.me = - self.mur * np.sqrt(self.energy)

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

    def fromstring(self, string):
        l = string.split()
        self.energy = float(l.pop(0))
        self.mur = np.array([float(l.pop(0)) for i in range(3)])
        self.muv = np.array([float(l.pop(0)) for i in range(3)])


class FakeExcitedStates(ExcitationList):
    def __init__(self, calculator):
        ExcitationList.__init__(self, calculator)

        energies = np.arange(1, 4)
        murs = np.eye(3) / energies
        muvs = np.eye(3)
        for energy, mur, muv in zip(energies, murs, muvs):
            self.append(FakeExcitation(energy, mur, muv))

        self.calculator = calculator
    
    def read(self, filename):
        """Read myself from a file"""
        with open(filename, 'r') as f:
            self.filename = filename
            n = int(f.readline().split()[0])
            for i in range(n):
                self.append(FakeExcitation(string=f.readline()))

    def write(self, fname):
        with open(fname, 'w') as f:
            print(len(self), file=f)
            for ex in self:
                 f.write(ex.outstring())
                 

lj = LennardJones(sigma=2**(1/6))
atoms = Atoms('H2', positions=[[0, 0 ,0], [0, 0, 1]])
atoms.set_calculator(lj)

print('??? 1a')

pz = Placzek(atoms, FakeExcitedStates, txt='-')
pz.run()
print('??? 2')

om = 0.1
pzi = pz.absolute_intensity(omega=om)[-1]
print(pzi)
#print(pz.exm_r)
