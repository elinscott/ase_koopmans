from ase.calculators.calculator import Calculator, all_changes
from ase.calculators.calculator import PropertyNotImplementedError


class MixingCalculator(Calculator):
    """Special calculator for combining multiple calculators.
    """

    def __init__(self, *calculators, **kwargs):
        """Implementation of mixing of calculators.

        calculators: list
            List of any arbitrary Calculators
        kwargs: dict
            Additional key value pairs for the `Calculator` class.
        """
        # TODO: initialisation should be happen as multiple arguments or as single list as argument?

        super().__init__(**kwargs)

        if kwargs.get('restart', None) or kwargs.get('ignore_bad_restart_file', None):
            raise NotImplementedError('The restart and the ignore_bad_restart_file options '
                                      'doesn\'t make any sense here.')

        if kwargs.get('label', None) or kwargs.get('directory', None):
            raise NotImplementedError('Please use the label and the directory options for the individual calculators.')

        for calc in calculators:
            if not isinstance(calc, Calculator):
                raise ValueError('All the calculators should be inherited form the ase\'s Calculator class')

        common_properties = set.intersection(*(set(calc.implemented_properties) for calc in calculators))
        self.implemented_properties = list(common_properties)

        if not self.implemented_properties:
            raise PropertyNotImplementedError('There are no common property implemented for the potentials!')

        self.calculators = calculators

    def calculate(self, atoms=None, properties=['energy'], system_changes=all_changes):
        """ Calculates all the specific property for each calculator and returns with the summed value.
        """

        super().calculate(atoms, properties, system_changes)

        if not set(properties).issubset(self.implemented_properties):
            raise PropertyNotImplementedError('Some of the requested property is not in the '
                                              'list of supported properties ({})'.format(self.implemented_properties))

    def reset(self):
        """Clear all previous results recursively from all fo the calculators."""
        super().reset()

        for calc in self.calculators:
            calc.reset()

    def __str__(self):
        # calculators = ', '.join(str(calc) for calc in self.calculators)
        calculators = ', '.join(calc.__class__.__name__ for calc in self.calculators)
        return '{}({})'.format(self.__class__.__name__, calculators)


class SumCalculator(MixingCalculator):
    """SumCalculator for combining multiple calculators.

    This calculator can be used when there are different calculators for the different chemical environment or
    for example during delta leaning. It works with a list of arbitrary calculators and evaluates them in sequence
    when it is required.
    The supported properties are the intersection of the implemented properties in each calculator.
    """

    def calculate(self, atoms=None, properties=['energy'], system_changes=all_changes):
        """ Calculates all the specific property for each calculator and returns with the summed value.
        """

        super().calculate(atoms, properties, system_changes)

        for calc in self.calculators:
            calc.calculate(atoms, properties, system_changes)

            for k in properties:
                if k not in self.results:
                    self.results[k] = calc.results[k]
                else:
                    self.results[k] += calc.results[k]


class LinearCombinationCalculator(MixingCalculator):
    """LinearCombinationCalculator for weighted summation of multiple calculators (for thermodynamic purposes)..
    """

    def __init__(self, *calculators, weights=None, **kwargs):
        """Implementation of sum of calculators.

        calculators: list
            List of any arbitrary Calculators
        weights: list of float
            List of any arbitrary Calculators
        kwargs: dict
            Additional key value pairs for the `Calculator` class.
        """

        super().__init__(*calculators, **kwargs)

        if weights is None:
            n = len(calculators)
            weights = [1 / n] * n  # average by default

        if len(weights) != len(calculators):
            raise ValueError('The length of the weights must be the same as the number of calculators!')

        self.weights = weights

    def calculate(self, atoms=None, properties=['energy'], system_changes=all_changes):
        """ Calculates all the specific property for each calculator and returns with the summed value.
        """

        super().calculate(atoms, properties, system_changes)

        for w, calc in zip(self.weights, self.calculators):
            calc.calculate(atoms, properties, system_changes)

            for k in properties:
                if k not in self.results:
                    self.results[k] = w * calc.results[k]
                else:
                    self.results[k] += w * calc.results[k]
