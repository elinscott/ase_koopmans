"""
Exceptions for the general error types that can occur either while
setting up the calculator, which requires constructing KIM API C++
objects, or while running a simulation
"""
from ase.calculators.calculator import CalculatorError


class KIMCalculatorError(CalculatorError):
    """
    Indicates an error occurred in initializing the calculator, e.g. due
    to incompatible combinations of argument values
    """

    pass


class KIMModelNotFound(CalculatorError):
    """
    Requested model cannot be found in any of the KIM API model
    collections
    """

    pass


class KIMInitializationError(CalculatorError):
    """
    KIM API Model or ComputeArguments objects could not be successfully
    created
    """

    pass


class KimpyError(CalculatorError):
    """
    A call to a kimpy function returned a non-zero error code
    """

    pass
