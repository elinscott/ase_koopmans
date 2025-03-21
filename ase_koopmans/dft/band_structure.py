from ase_koopmans.spectrum.band_structure import *  # noqa: F401,F403
import warnings
from numpy import VisibleDeprecationWarning

warnings.warn("ase_koopmans.dft.band_structure has been moved to "
              "ase_koopmans.spectrum.band_structure. Please_koopmans update your "
              "scripts; this alias will be removed in a future "
              "version of ASE.",
              VisibleDeprecationWarning)
