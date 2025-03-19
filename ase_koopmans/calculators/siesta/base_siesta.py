# File solely for gently deprecating this Base_koopmansSiesta class.
import warnings
import numpy as np
from ase_koopmans.calculators.siesta.siesta import Siesta


class Base_koopmansSiesta(Siesta):
    def __init__(self, *args, **kwargs):
        warnings.warn(
            "The Base_koopmansSiesta calculator class will no longer be supported. "
            "Use 'ase_koopmans.calculators.siesta.Siesta in stead.",
            np.VisibleDeprecationWarning)
        Siesta.__init__(self, *args, **kwargs)
