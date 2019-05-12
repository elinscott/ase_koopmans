# Importing kimpy here is solely done to force ASE's CI to skip
# this calculator
import kimpy as _kimpy
#from ase.calculators.kim.kim import KIM
from .kim import KIM

# Use the kimpy module imported to avoid flake8 warning (currently
# there's no actual version requirements on kimpy)
_ = _kimpy

__all__ = ["KIM"]
