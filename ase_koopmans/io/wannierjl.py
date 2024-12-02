"""
This module defines I/O routines with WannierJL files

In reality, WannierJL does not rely on input and output files
"""

from ase_koopmans import Atoms
from ase_koopmans.calculators.wannierjl import WannierJL


def read_wannierjl_out(fd):
    """
    Dummy function that pretends to read WannierJL output files

    Parameters
    ----------
    fd : file|str
        A file like object or filename

    Yields
    ------
    structure : atoms
        An Atoms object with an attached WannierJL Calculator containing
        any parsed results
    """

    structure = Atoms()
    calc = WannierJL(atoms=structure)
    structure.calc = calc

    yield structure