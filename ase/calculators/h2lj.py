import numpy as np
from ase import Atoms
from ase.units import invcm
from ase.calculators.calculator import all_changes
from ase.calculators.lj import LennardJones
from ase.calculators.excitation import ExcitationList, Excitation

# data from:
# https://webbook.nist.gov/cgi/cbook.cgi?ID=C1333740&Mask=1000#Diatomic
#     X        B        C
Re = [0.74144, 1.2928,  1.0327,  1.0327 ]  # eq. bond length
ome= [4401.21, 1358.09, 2443.77, 2443.77]  # vibrational frequency
ome = np.array(ome) * invcm
Ee = [0,       91700.0, 100089.9, 100089.9]  # transition energy at Re[0]
Ee = np.array(Ee) * invcm

# LJ parameters
sigma = np.array(Re) * 2**(-1 / 6)
k = ome**2 * 0.5 * (4401.21 / 226.2)**2  # XXXX correct factor?
epsilon = k * sigma**2 / 72 / 2**(1 / 3)
Ee -= epsilon[0]

def H2LJ():
    """Return Lennard-Jones H2 with calculator attached."""
    atoms = Atoms('H2', positions=np.zeros((2, 3)))
    atoms[1].position[2] = Re[0]
    atoms.set_calculator(H2ljState(0))
    atoms.get_potential_energy()
    return atoms

class H2ljState(LennardJones):
    """H2 ground state as Lennard-Jones potential"""
    def __init__(self, index):
        LennardJones.__init__(self,
                              sigma=sigma[index], epsilon=epsilon[index])

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=all_changes):
        if atoms is not None:
            assert len(atoms) == 2
        LennardJones.calculate(self, atoms, properties, system_changes)


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


class H2LJExcitedStates(ExcitationList):
    """First singlet excited state of H2 as Lennard-Jones potentials"""
    def __init__(self, calculator):
        ExcitationList.__init__(self, calculator)

        h2lj = H2LJ()

        # molecular axis
        atoms = calculator.get_atoms()
        vr = atoms[1].position - atoms[0].position
        r = np.linalg.norm(vr)
        hr = vr / r
        # perpendicular axes
        vrand = np.random.rand(3)
        hx = np.cross(hr, vrand)
        hx /= np.linalg.norm(hx)
        hy = np.cross(hr, hx)
        hy /= np.linalg.norm(hy)

        # central me value and rise
        hvec = [None, hr, hx, hy]
        mc = [0, 0.9, 0.8, 0.8]
        mr = [0, 1.0, 0.5, 0.5]
            
        for i in range(1, 4):
            calc = H2ljState(i)
            calc.calculate(calculator.get_atoms())
            energy = Ee[i] + calc.get_potential_energy()
            calc.calculate(h2lj)
            energy -= calc.get_potential_energy()
            
            mur = hvec[i] * (mc[i] + (r - Re[0]) * mr[i])
            muv = mur 

            self.append(FakeExcitation(energy, mur, muv))

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
