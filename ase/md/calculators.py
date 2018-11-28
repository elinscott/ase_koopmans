import numpy as np
from ase import units
from ase.calculators.calculator import Calculator, all_changes


class MixedCalculator(Calculator):
    """
    Mixing of two calculators with different weights

    H = weight1 * H1 + weight2 * H2

    Parameters
    ----------
    calc1 : ASE-calculator
    calc2 : ASE-calculator
    weight1 : float
        weight for calculator 1
    weight2 : float
        weight for calculator 2
    """

    implemented_properties = ['forces', 'energy']

    def __init__(self, calc1, calc2, weight1, weight2):
        Calculator.__init__(self)
        self.calc1 = calc1
        self.calc2 = calc2
        self.weight1 = weight1
        self.weight2 = weight2

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=all_changes):
        Calculator.calculate(self, atoms, properties, system_changes)

        atoms_tmp = atoms.copy()
        atoms_tmp.set_calculator(self.calc1)
        forces1 = atoms_tmp.get_forces()
        energy1 = atoms_tmp.get_potential_energy()

        atoms_tmp = atoms.copy()
        atoms_tmp.set_calculator(self.calc2)
        forces2 = atoms_tmp.get_forces()
        energy2 = atoms_tmp.get_potential_energy()

        self.results['forces'] = self.weight1*forces1 + self.weight2*forces2
        self.results['energy'] = self.weight1*energy1 + self.weight2*energy2
        self.results['energy_contributions'] = (energy1, energy2)

    def get_energy_contributions(self):
        """ Return the potential energy from calc1 and calc2 respectively """
        return self.results['energy_contributions']


class EinsteinCalculator(Calculator):
    """
    Einstein crystal calculator

    Parameters
    ----------
    ideal_positions : array
        array of the ideal crystal positions
    k : float
        spring constant
    """
    implemented_properties = ['forces', 'energy']

    def __init__(self, ideal_positions, k):
        Calculator.__init__(self)
        self.ideal_positions = ideal_positions.copy()
        self.k = k

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=all_changes):
        Calculator.calculate(self, atoms, properties, system_changes)
        energy, forces = self.compute_energy_and_forces(atoms)
        self.results['energy'], self.results['forces'] = energy, forces

    def compute_energy_and_forces(self, atoms):
        disps = atoms.positions - self.ideal_positions
        forces = - self.k * disps
        energy = sum(self.k / 2.0 * np.linalg.norm(disps, axis=1)**2)
        return energy, forces


def get_einstein_free_energy(k, m, T, method='classical'):
    """ Get free energy (per atom) for an Einstein crystal.

    Free energy of a Einstein solid given by classical (1) or QM (2)
    1.    F_E = 3NkbT log( hw/kbT )
    2.    F_E = 3NkbT log( 1-exp(hw/kbT) ) + zeropoint

    Parameters
    -----------
    k : float
        spring constant (eV/A^2)
    m : float
        mass (grams/mole or AMU)
    T : float
        temperature (K)
    method : str
        method for free energy computation, classical or QM.

    Returns
    --------
    float
        free energy of the Einstein crystal (eV/atom)
    """
    assert method in ['classical', 'QM']

    hbar = units._hbar * units.J  # eV/s
    m = m / units.kg              # mass kg
    k = k * units.m**2 / units.J  # spring constant J/m2
    omega = np.sqrt(k / m)        # angular frequency 1/s

    if method == 'classical':
        F_einstein = 3 * units.kB * T * np.log(hbar*omega/(units.kB*T))
    elif method == 'QM':
        F_einstein = 3 * units.kB * T * np.log(1.0 - np.exp(-hbar*omega/(units.kB*T))) + 1.5 * hbar*omega

    return F_einstein
