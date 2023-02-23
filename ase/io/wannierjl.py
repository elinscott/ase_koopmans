"""
This module defines I/O routines with WannierJL files.
Note that WannierJL uses Wannier90 input files.
"""

from ase import Atoms
from ase.utils import basestring
from ase.calculators.wannierjl import WannierJL


def read_wannierjl_out(fd):
    """
    Reads wannierjl output files

    Parameters
    ----------
    fd : file|str
        A file like object or filename

    Yields
    ------
    structure : atoms
        An Atoms object with an attached SinglePointCalculator containing
        any parsed results
    """

    if isinstance(fd, basestring):
        fd = open(fd, 'rU')

    flines = fd.readlines()

    structure = Atoms()

    for i, line in enumerate(flines):
        pass

    calc = WannierJL(atoms=structure)

    structure.calc = calc

    yield structure