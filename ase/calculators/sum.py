import numpy as np

from ase.calculators.calculator import Calculator, all_properties, all_changes
from ase.calculators.calculator import PropertyNotImplementedError


class SumCalculator(Calculator):
    """Special calculator for combining multiple calculators.

    This calculator can be used when there are different calculators for the different chemical environment or
    for example during delta leaning. It works with a list of arbitrary calculators and evaluates them in sequence
    when it is required.
    The supported properties are the intersection ot the implemented properties in each calculator.
    """

    def __init__(self, *calculators, **kwargs):
        """Implementation of sum of calculators.

        calculators: list
            List of any arbitrary Calculators
        kwargs: dict
            Additional key value pairs for th `Calculator` class.
        """
        # TODO: initialisation should be happen as multiple arguments or as single list as argument?

        super(SumCalculator, self).__init__(**kwargs)

        if kwargs.get('restart', None) or kwargs.get('ignore_bad_restart_file', None):
            raise NotImplementedError(
                'The restart and the ignore_bad_restart_file options doesn\'t make any sense here.')

        if kwargs.get('label', None) or kwargs.get('directory', None):
            raise NotImplementedError('Please use the label and the directory options for the individual calculators.')

        for calc in calculators:
            if not isinstance(calc, Calculator):
                raise ValueError('All the calculators should be inherited form the ase\'s Calculator class')

        common_properties = set.intersection(*(set(calc.implemented_properties) for calc in calculators))
        self.implemented_properties = list(common_properties)

        if not self.implemented_properties:
            raise NotImplementedError('There are no common property implemented for the potentials!')

        self.calculators = calculators

    def calculate(self, atoms=None, properties=['energy'], force_consistent=False, system_changes=all_changes):
        """ Calculates all the specific property for each calculator and returns with the summed value.
        """

        super(SumCalculator, self).calculate(atoms, properties, system_changes)

        if not set(properties).issubset(self.implemented_properties):
            raise NotImplementedError(
                'Some ot the requested property is not in the list of supported properties ({})'.format(
                    self.implemented_properties))

        for calc in self.calculators:
            calc.calculate(atoms, properties, system_changes)

            for k in properties:
                if k not in self.results:
                    self.results[k] = calc.results[k]
                else:
                    self.results[k] += calc.results[k]

    def reset(self):
        """Clear all previous results recursively from all fo the calculators."""
        super(SumCalculator, self).reset()

        for calc in self.calculators:
            calc.reset()

    def __str__(self):
        # calculators = ', '.join(str(calc) for calc in self.calculators)
        calculators = ', '.join(calc.__class__.__name__ for calc in self.calculators)
        return '{}({})'.format(self.__class__.__name__, calculators)


if __name__ == '__main__':
    from ase.build import fcc111
    from ase.calculators.emt import EMT
    from ase.constraints import FixAtoms

    # Calculate reference values:
    atoms = fcc111('Cu', (2, 2, 1), vacuum=10.)
    atoms[0].x += 0.2

    # First run the test with EMT similarly to the test of the single point calculator.
    calc = EMT()
    atoms.set_calculator(calc)
    forces = atoms.get_forces()

    # Alternative ways to associate a calaculator with an atoms object.
    atoms1 = atoms.copy()
    calc1 = SumCalculator(EMT(), EMT())
    atoms1.set_calculator(calc1)

    atoms2 = atoms.copy()
    calc2 = SumCalculator(EMT(), EMT(), atoms=atoms2)

    # Check the results.
    assert np.isclose(2 * forces, atoms1.get_forces()).all()
    assert np.isclose(2 * forces, atoms2.get_forces()).all()

    # testing  step
    atoms1[0].x += 0.2
    assert not np.isclose(2 * forces, atoms1.get_forces()).all()

    # Check constraints
    atoms1.set_constraint(FixAtoms(indices=[atom.index for atom in atoms]))
    assert np.isclose(0, atoms1.get_forces()).all()

    # Representation
    print(calc1)
