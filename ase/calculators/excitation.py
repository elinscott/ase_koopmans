"""Excitation lists base classes

"""
import numpy as np

from ase.parallel import world
from ase.utils import convert_string_to_fd
from ase.units import Hartree

class ExcitationList(list):

    """General Excitation List class.

    """

    def __init__(self, calculator=None, txt='-'):
        """
        Parameters
        ----------
        calculator: object or string
          if calculator is a string: read 
          else: calculate
        """
        # initialise empty list
        list.__init__(self)

        if not txt and hasattr(calculator, 'log'):
            txt = calculator.log.fd
        self.txt = convert_string_to_fd(txt, world)

        if isinstance(calculator, str):
            # initialize from a file
            self.calculator = None
            self.read(calculator)
        else:
            self.calculator = calculator
            self.calculate()

    def get_calculator(self):
        return self.calculator

    def get_energies(self):
        """Get excitation energies in Hartrees"""
        el = []
        for ex in self:
            el.append(ex.get_energy())
        return np.array(el)

    def get_trk(self):
        """Evaluate the Thomas Reiche Kuhn sum rule"""
        trkm = np.zeros((3))
        for ex in self:
            me = ex.get_dipole_me()
            trkm += ex.get_energy() * (me.real ** 2 + me.imag ** 2)
        return 2. * trkm  # scale to get the number of electrons XXX spinpol ?

    def get_polarizabilities(self, lmax=7):
        """Calculate the Polarisabilities
        see Jamorski et al. J. Chem. Phys. 104 (1996) 5134"""
        S = np.zeros((lmax + 1))
        for ex in self:
            e = ex.get_energy()
            f = ex.get_oscillator_strength()[0]
            for l in range(lmax + 1):
                S[l] += e ** (-2 * l) * f
        return S

    def set_calculator(self, calculator):
        self.calculator = calculator

    def __truediv__(self, x):
        return self.__mul__(1. / x)

    __div__ = __truediv__

    def __rmul__(self, x):
        return self.__mul__(x)

    def __mul__(self, x):
        """Multiply with a number"""
        if isinstance(x, (float, int)):
            result = self.__class__()
            result.dtype = self.dtype
            for kss in self:
                result.append(x * kss)
            return result
        else:
            return RuntimeError('not a number')

    def __sub__(self, other):
        result = self.__class__()
        result.dtype = self.dtype
        assert(len(self) == len(other))
        for kss, ksso in zip(self, other):
            result.append(kss - ksso)
        return result

    def __str__(self):
        string = '# ' + str(type(self))
        if len(self) != 0:
            string += ', %d excitations:' % len(self)
        string += '\n'
        for ex in self:
            string += '#  ' + ex.__str__() + '\n'
        return string


class Excitation:
    """Base class for a single excited state"""
    def get_energy(self):
        """Get the excitations energy relative to the ground state energy
        in Hartrees.
        """
        return self.energy

    def get_dipole_me(self, form='r'):
        """return the excitations dipole matrix element
        including the occupation factor sqrt(fij)"""
        if form == 'r':
            # length form
            return self.me / np.sqrt(self.energy)
        elif form == 'v':
            # velocity form
            return - np.sqrt(self.fij) * self.muv
        else:
            raise RuntimeError('Unknown form >' + form + '<')

    def get_dipole_tensor(self, form='r'):
        """Return the "oscillator strength tensor"

        self.me is assumed to be::

          form='r': sqrt(f * E) * <I|r|J>,
          form='v': sqrt(f / E) * <I|d/(dr)|J>

        for f = multiplicity, E = transition energy and initial and
        final states::

          |I>, |J>
        """

        if form == 'r':
            # length form
            me = self.me
        elif form == 'v':
            # velocity form
            me = self.muv * np.sqrt(self.fij * self.energy)
        else:
            raise RuntimeError('Unknown form >' + form + '<')
        
        return 2 * np.outer(me, me.conj())

    def get_oscillator_strength(self, form='r'):
        """Return the excitations dipole oscillator strength."""
        me2_c = self.get_dipole_tensor(form).diagonal().real
        return np.array([np.sum(me2_c) / 3.] + me2_c.tolist())

    def get_rotatory_strength(self, form='r', units='cgs'):
        """Return rotatory strength"""
        if self.magn is None:
            raise RuntimeError('Magnetic moment not available.')

        if units == 'cgs':
            # 10^-40 esu cm erg / G
            # = 3.33564095 * 10^-15 A^2 m^3 s
            # conversion factor after
            # T. B. Pedersen and A. E. Hansen,
            # Chem. Phys. Lett. 246 (1995) 1
            # pre = 471.43
            # From TurboMole
            pre = 64604.8164
        elif units == 'a.u.':
            pre = 1.
        else:
            raise RuntimeError('Unknown units >' + units + '<')

        if form == 'r':
            # length form
            mu = self.mur
        elif form == 'v':
            # velocity form
            mu = self.muv
        else:
            raise RuntimeError('Unknown form >' + form + '<')

        return pre * np.dot(mu, self.magn)

    def __str__(self):
        string = '# <' + self.__class__.__name__ + '> '
        string += 'eji=%g[eV]' % self.energy
        if self.me.dtype == float:
            string += ' (%g,%g,%g)' % (self.me[0], self.me[1], self.me[2])
        else:
            string += ' kpt={0:d} w={1:g}'.format(self.k, self.weight)
            string += ' ('
            # use velocity form
            s = - np.sqrt(self.energy * self.fij)
            for c, m in enumerate(s * self.me):
                string += '{0.real:.5e}{0.imag:+.5e}j'.format(m)
                if c < 2:
                    string += ','
            string += ')'
        return string

    def set_energy(self, E):
        """Set the excitations energy relative to the ground state energy"""
        self.energy = E
