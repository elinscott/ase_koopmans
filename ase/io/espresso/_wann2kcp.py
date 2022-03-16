"""Reads wann2kcp files.

"""

from ase.calculators.espresso import Wann2KCP
from ._x2y import read_x2y_in, write_x2y_in, read_x2y_out


def read_wann2kcp_in(fileobj):
    # Parse a wann2kcp input file, 'w2kcp', '.wki'
    return read_x2y_in(fileobj, Wann2KCP)


def write_wann2kcp_in(fd, atoms, **kwargs):
    # Create an input file for wann2kcp.
    return write_x2y_in(fd, atoms, **kwargs)


def read_wann2kcp_out(fd):
    # Reads wann2kcp output files
    yield from read_x2y_out(fd, Wann2KCP)
