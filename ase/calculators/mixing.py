from ase.calculators.calculator import Calculator, all_changes
from ase.calculators.calculator import PropertyNotImplementedError


class LinearCombinationCalculator(Calculator):
    """LinearCombinationCalculator for weighted summation of multiple calculators (for thermodynamic purposes)..
    """

    def __init__(self, calcs=None, weights=None, atoms=None):
        """Implementation of sum of calculators.

        calcs: list
            List of any arbitrary Calculators
        weights: list of float
            List of any arbitrary Calculators
        """

        super().__init__(atoms=atoms)

        for calc in calcs:
            if not isinstance(calc, Calculator):
                raise ValueError('All the calculators should be inherited form the ase\'s Calculator class')

        common_properties = set.intersection(*(set(calc.implemented_properties) for calc in calcs))
        self.implemented_properties = list(common_properties)

        if not self.implemented_properties:
            raise PropertyNotImplementedError('There are no common property implemented for the potentials!')

        if len(weights) != len(calcs):
            raise ValueError('The length of the weights must be the same as the number of calculators!')

        self.calcs = calcs
        self.weights = weights

    def calculate(self, atoms=None, properties=['energy'], system_changes=all_changes):
        """ Calculates all the specific property for each calculator and returns with the summed value.
        """

        super().calculate(atoms, properties, system_changes)

        if not set(properties).issubset(self.implemented_properties):
            raise PropertyNotImplementedError('Some of the requested property is not in the '
                                              'list of supported properties ({})'.format(self.implemented_properties))

        for w, calc in zip(self.weights, self.calcs):
            calc.calculate(atoms, properties, system_changes)

            for k in properties:
                if k not in self.results:
                    self.results[k] = w * calc.results[k]
                else:
                    self.results[k] += w * calc.results[k]

    def reset(self):
        """Clear all previous results recursively from all fo the calculators."""
        super().reset()

        for calc in self.calcs:
            calc.reset()

    def __str__(self):
        # calculators = ', '.join(str(calc) for calc in self.calculators)
        calculators = ', '.join(calc.__class__.__name__ for calc in self.calcs)
        return '{}({})'.format(self.__class__.__name__, calculators)


class SumCalculator(LinearCombinationCalculator):
    """SumCalculator for combining multiple calculators.

    This calculator can be used when there are different calculators for the different chemical environment or
    for example during delta leaning. It works with a list of arbitrary calculators and evaluates them in sequence
    when it is required.
    The supported properties are the intersection of the implemented properties in each calculator.
    """

    def __init__(self, calcs=None, atoms=None):
        """Implementation of sum of calculators.

        calcs: list
            List of any arbitrary Calculators
        """

        weights = [1.] * len(calcs)
        super().__init__(calcs=calcs, weights=weights, atoms=atoms)


class AverageCalculator(LinearCombinationCalculator):
    """AverageCalculator for weighted summation of multiple calculators (for thermodynamic purposes)..
    """

    def __init__(self, calcs=None, atoms=None):
        """Implementation of average of calculators.

        calcs: list
            List of any arbitrary Calculators
        """

        n = len(calcs)
        weights = [1 / n] * n
        super().__init__(calcs=calcs, weights=weights, atoms=atoms)
