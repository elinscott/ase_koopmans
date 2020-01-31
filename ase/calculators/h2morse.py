import numpy as np
from ase import Atoms
from ase.units import invcm, Ha
from ase.data import atomic_masses
from ase.calculators.calculator import all_changes
from ase.calculators.morse import MorsePotential
from ase.calculators.excitations import ExcitationList, Excitation


# data from:
# https://webbook.nist.gov/cgi/cbook.cgi?ID=C1333740&Mask=1000#Diatomic
#     X        B       C       C
Re = [0.74144, 1.2928, 1.0327, 1.0327]  # eq. bond length
ome = [4401.21, 1358.09, 2443.77, 2443.77]  # vibrational frequency
ome = np.array(ome)
# electronic transition energy
Etrans = [0, 91700.0, 100089.9, 100089.9]
Etrans = np.array(Etrans) * invcm

# dissociation energy
# GS: https://aip.scitation.org/doi/10.1063/1.3120443
De = np.ones(4) * 36118.069 * invcm
# B, C separated energy E(1s) - E(2p)
De[1:] += Ha / 2 - Ha / 8
De -= Etrans

# Morse parameter
m = atomic_masses[1] * 0.5  # reduced mass
# XXX find scaling factor
rho0 = Re * ome * invcm * np.sqrt(m / 2 / De) * 4401.21 / 284.55677429605862


def H2Morse(state=0):
    """Return H2 as a Morse-Potential with calculator attached."""
    atoms = Atoms('H2', positions=np.zeros((2, 3)))
    atoms[1].position[2] = Re[state]
    atoms.set_calculator(H2MorseState(state))
    atoms.get_potential_energy()
    return atoms


class H2MorseState(MorsePotential):
    """H2 ground or excited state as Morse potential"""
    def __init__(self, state):
        MorsePotential.__init__(self,
                                epsilon=De[state],
                                r0=Re[state], rho0=rho0[state])

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=all_changes):
        if atoms is not None:
            assert len(atoms) == 2
        MorsePotential.calculate(self, atoms, properties, system_changes)

        # determine 'wave functions' including
        # Berry phase (arbitrary sign) and
        # random orientation of wave functions perpendicular
        # to the molecular axis
        
        # molecular axis
        atoms = self.calculator.get_atoms()
        vr = atoms[1].position - atoms[0].position
        r = np.linalg.norm(vr)
        hr = vr / r
        # perpendicular axes
        vrand = np.random.rand(3)
        hx = np.cross(hr, vrand)
        hx /= np.linalg.norm(hx)
        hy = np.cross(hr, hx)
        hy /= np.linalg.norm(hy)
        wfs = [1, hr, hx, hy]
        # Berry phase
        berry = (-1)**np.random.randint(0, 2, 4)
        self.wfs = [wf * b for wf, b in zip(wfs, berry)]
        
    def write(filename, option=None):
        """write calculated state to a file"""
        with open(filename, 'w') as f:
            for wf in self.wfs:
                f.write(str(wf) + '\n')
 
class H2MorseExcitedStates(ExcitationList):
    """First singlet excited state of H2 as Lennard-Jones potentials"""
    def __init__(self, calculator):
        ExcitationList.__init__(self, calculator)

    def calculate(self):
        """Caclculate excitation spectrum"""
        # molecular axis
        atoms = self.calculator.get_atoms()
        vr = atoms[1].position - atoms[0].position
        r = np.linalg.norm(vr)
        hr = vr / r
        # perpendicular axes
        vrand = np.random.rand(3)
        hx = np.cross(hr, vrand)
        hx /= np.linalg.norm(hx)
        hy = np.cross(hr, hx)
        hy /= np.linalg.norm(hy)

        # matrix element directions
        hvec = [None, hr, hx, hy]
        # central me value and rise, unit Bohr
        # from DOI: 10.1021/acs.jctc.9b00584
        mc = [0, 0.8, 0.7, 0.7]
        mr = [0, 1.0, 0.5, 0.5]

        calc = H2MorseState(0)
        calc.calculate(atoms)
        E0 = calc.get_potential_energy()
        for i in range(1, 4):
            energy = Ha * (0.5 - 1. / 8) - E0
            calc = H2MorseState(i)
            calc.calculate(atoms)
            energy += calc.get_potential_energy()

            mur = hvec[i] * (mc[i] + (r - Re[0]) * mr[i])
            muv = mur

            self.append(BasicExcitation(energy, i, mur, muv))

    def read(self, filename):
        """Read myself from a file"""
        with open(filename, 'r') as f:
            self.filename = filename
            n = int(f.readline().split()[0])
            for i in range(n):
                self.append(BasicExcitation(string=f.readline()))

    def write(self, fname):
        with open(fname, 'w') as f:
            f.write('{0}\n'.format(len(self)))
            for ex in self:
                f.write(ex.outstring())


class BasicExcitation(Excitation):
    def __init__(self, energy=None, index=None,
                 mur=None, muv=None, magn=None, string=None):
        if string is not None:
            self.fromstring(string)
        else:
            self.energy = energy
            self.index = index
            self.mur = mur
            self.muv = muv
            self.magn = magn
        self.fij = 1.
        self.me = - self.mur * np.sqrt(self.energy / Ha)

    def __eq__(self, other):
        """Considered to be equal when their indices are equal."""
        return self.index == other.index

    def __hash__(self):
        """Hash similar to __eq__"""
        if not hasattr(self, 'hash'):
            self.hash = hash(self.index)
        return self.hash

    def outstring(self):
        string = '{0:g}  {1}  '.format(self.energy, self.index)

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
        energy = float(l.pop(0))
        index = int(l.pop(0))
        mur = np.array([float(l.pop(0)) for i in range(3)])
        muv = np.array([float(l.pop(0)) for i in range(3)])
        self.__init__(energy, index, mur, muv)
