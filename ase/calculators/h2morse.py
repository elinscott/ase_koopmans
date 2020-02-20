import numpy as np
from ase import Atoms
from ase.units import invcm, Ha
from ase.data import atomic_masses
from ase.calculators.calculator import all_changes
from ase.calculators.morse import MorsePotential
from ase.calculators.excitation_list import Excitation, ExcitationList

"""The H2 molecule represented by Morse-Potentials for
gound and first 3 excited singlet states B + C(doubly degenerate)"""


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
    atoms.calc = H2MorseState(state)
    atoms.get_potential_energy()
    return atoms


class H2MorseState(MorsePotential):
    """H2 ground or excited state as Morse potential"""
    def __init__(self, state):
        if isinstance(state, str):
            self.read(state)
        else:
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

    def read(self, filename):
        with open(filename) as f:
            self.wfs = [int(f.readline().split()[0])]
            for i in range(1, 4):
                self.wfs.append(
                    np.array([float(x)
                              for x in f.readline().split()[:4]]))
        
    def write(self, filename, option=None):
        """write calculated state to a file"""
        with open(filename, 'w') as f:
            f.write('{}\n'.format(self.wfs[0]))
            for wf in self.wfs[1:]:
                f.write('{0:g} {1:g} {2:g}\n'.format(*wf))

    def overlap(self, other):
        ov = np.zeros((4, 4))
        ov[0, 0] = self.wfs[0] * other.wfs[0]
        wfs = np.array(self.wfs[1:])
        owfs = np.array(other.wfs[1:])
        ov[1:, 1:] = np.dot(wfs, owfs.T)
        return ov


class H2MorseExcitedStatesCalculator():
    """First singlet excited state of H2 as Lennard-Jones potentials"""
    def __init__(self, gscalculator=None, nstates=3, txt='-'):
        """
        Parameters
        ----------
        gscalculator: object
          Calculator for ground state energies
        nstates: int
          Numer of states to calculate 0 < nstates < 4, default 3
        txt:
          output channel, default '-'
        """
        assert nstates > 0 and nstates < 4
        self.nstates = nstates
        self.gscalc = gscalculator

    def calculate(self, atoms=None):
        """Calculate excitation spectrum

        Parameters
        ----------
        atoms: Ase atoms object
           Default None
        """
        # central me value and rise, unit Bohr
        # from DOI: 10.1021/acs.jctc.9b00584
        mc = [0, 0.8, 0.7, 0.7]
        mr = [0, 1.0, 0.5, 0.5]

        if atoms is None:
            atoms = self.gscalc.atoms
        
        if self.gscalc is not None:
            cgs = self.gscalc
        else:
            cgs = atoms.calc
        r = atoms.get_distance(0, 1)
        E0 = cgs.get_potential_energy(atoms)
        
        exl = H2MorseExcitedStates()
        for i in range(1, self.nstates + 1):
            hvec = cgs.wfs[0] * cgs.wfs[i]
            energy = Ha * (0.5 - 1. / 8) - E0
            calc = H2MorseState(i)
            calc.calculate(atoms)
            energy += calc.get_potential_energy()

            mur = hvec * (mc[i] + (r - Re[0]) * mr[i])
            muv = mur

            exl.append(H2Excitation(energy, i, mur, muv))
        return exl


class H2MorseExcitedStates(ExcitationList):
    """First singlet excited state of H2 as Lennard-Jones potentials"""
    def __init__(self, filename=None, nstates=3):
        self.nstates = nstates
        ExcitationList.__init__(self, filename)

    def overlap(self, ov_nn, other):
        return (ov_nn[1:len(self) + 1, 1:len(self) + 1] *
                ov_nn[0, 0])

    def read(self, filename):
        """Read myself from a file"""
        with open(filename, 'r') as f:
            self.filename = filename
            n = int(f.readline().split()[0])
            for i in range(min(n, self.nstates)):
                self.append(H2Excitation.fromstring(f.readline()))

    def write(self, fname):
        with open(fname, 'w') as f:
            f.write('{0}\n'.format(len(self)))
            for ex in self:
                f.write(ex.outstring())


class H2Excitation(Excitation):
    def __eq__(self, other):
        """Considered to be equal when their indices are equal."""
        return self.index == other.index

    def __hash__(self):
        """Hash similar to __eq__"""
        if not hasattr(self, 'hash'):
            self.hash = hash(self.index)
        return self.hash


class H2MorseExcitedStatesAndCalculator(
        H2MorseExcitedStatesCalculator, H2MorseExcitedStates):
    """Traditional joined object"""
    def __init__(self, calculator=None, nstates=3):
        if isinstance(calculator, str):
            H2MorseExcitedStates.__init__(self, calculator, nstates)
        else:
            excalc = H2MorseExcitedStatesCalculator(calculator, nstates)
            exlist = excalc.calculate()
            H2MorseExcitedStates.__init__(self, nstates=nstates)
            for ex in exlist:
                self.append(ex)
