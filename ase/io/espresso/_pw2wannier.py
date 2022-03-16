"""Reads pw2wannier files.

"""

from ase.calculators.espresso import PW2Wannier
from ._x2y import read_x2y_in, write_x2y_in, read_x2y_out


def read_pw2wannier_in(fileobj):
    # Parse a pw2wannier input file, 'pw2wan', '.p2wi'
    return read_x2y_in(fileobj, PW2Wannier)


def write_pw2wannier_in(fd, atoms, **kwargs):
    # Create an input file for pw2wannier.
    return write_x2y_in(fd, atoms, **kwargs)


def read_pw2wannier_out(fd):
    # Reads pw2wannier output files
    yield from read_x2y_out(fd, PW2Wannier)
