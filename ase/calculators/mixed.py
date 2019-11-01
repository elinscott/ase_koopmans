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

        for prop in properties:
            atoms_tmp = self.atoms.copy()
            prop1 = self.calc1.get_property(prop, atoms_tmp)
            atoms_tmp = self.atoms.copy()
            prop2 = self.calc2.get_property(prop, atoms_tmp)
            self.results[prop] = self.weight1 * prop1 + self.weight2 * prop2

            if prop == 'energy':
                self.results['energy_contributions'] = (prop1, prop2)

    def get_energy_contributions(self):
        """ Return the potential energy from calc1 and calc2 respectively """
        self.calculate(properties=['energy'])
        return self.results['energy_contributions']
