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
