import kimpy
from .kim import KIM

# Ensure minimal version requirement of kimpy is satsified
assert [int(x) for x in kimpy.__version__.split(".")] >= [0, 3, 2]

__all__ = ["KIM"]
