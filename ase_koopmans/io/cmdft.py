import numpy as np

from ase_koopmans.atom import Atom
from ase_koopmans.atoms import Atoms
from ase_koopmans.units import Bohr


def read_cmdft(fileobj):
    lines = fileobj.readlines()
    del lines[0]
    finished = False
    s = Atoms()
    while not finished:
        w = lines.pop(0).split()
        if w[0].startswith('"'):
            position = Bohr * np.array([float(w[3]), float(w[4]), float(w[5])])
            s.append(Atom(w[0].replace('"', ''), position))
        else:
            finished = True

    yield s
