import warnings
from collections import Counter

import numpy as np

from ase.data import atomic_numbers, chemical_symbols
from ase.utils.formula import Formula


def string2symbols(s):
    """Convert string to list of chemical symbols."""
    return list(Formula(s))


def symbols2numbers(symbols):
    if isinstance(symbols, str):
        symbols = string2symbols(symbols)
    numbers = []
    for s in symbols:
        if isinstance(s, str):
            numbers.append(atomic_numbers[s])
        else:
            numbers.append(int(s))
    return numbers


class Symbols:
    def __init__(self, numbers):
        self.numbers = numbers

    @classmethod
    def fromsymbols(cls, symbols):
        numbers = symbols2numbers(symbols)
        return cls(np.array(numbers))

    def formula(self):
        """Formula object."""
        return Formula.from_dict(Counter(self.get_chemical_symbols))

    def __getitem__(self, key):
        num = self.numbers[key]
        if np.isscalar(num):
            return chemical_symbols[num]
        return Symbols(num)

    def __setitem__(self, key, value):
        numbers = symbols2numbers(value)
        if len(numbers) == 1:
            numbers = numbers[0]
        self.numbers[key] = numbers

    def __len__(self):
        return len(self.numbers)

    def __str__(self):
        return self.get_chemical_formula('reduce')

    def __repr__(self):
        return 'Symbols(\'{}\')'.format(self)

    def __eq__(self, obj):
        if not hasattr(obj, '__len__'):
            return False

        try:
            symbols = Symbols.fromsymbols(obj)
        except Exception:
            # Typically this would happen if obj cannot be converged to
            # atomic numbers.
            return False
        return self.numbers == symbols.numbers

    def get_chemical_formula(self, mode='hill', empirical=False):
        """Get chemical formula.

        See documentation of ase.atoms.Atoms.get_chemical_formula()."""
        if mode in ('reduce', 'all') and empirical:
            warnings.warn("Empirical chemical formula not available "
                          "for mode '{}'".format(mode))

        if len(self) == 0:
            return ''

        numbers = self.numbers

        if mode == 'reduce':
            n = len(numbers)
            changes = np.concatenate(([0], np.arange(1, n)[numbers[1:] !=
                                                           numbers[:-1]]))
            symbols = [chemical_symbols[e] for e in numbers[changes]]
            counts = np.append(changes[1:], n) - changes

            tokens = []
            for s, c in zip(symbols, counts):
                tokens.append(s)
                if c > 1:
                    tokens.append(str(c))
            formula = ''.join(tokens)
        elif mode == 'all':
            formula = ''.join([chemical_symbols[n] for n in numbers])
        else:
            symbols = [chemical_symbols[Z] for Z in numbers]
            f = Formula('', [(symbols, 1)])
            if empirical:
                f, _ = f.reduce()
            if mode in {'hill', 'metal'}:
                formula = f.format(mode)
            else:
                raise ValueError(
                    "Use mode = 'all', 'reduce', 'hill' or 'metal'.")

        return formula
