import numpy as np
from ase.geometry.bravais import (recognize_canonical_cell,
                                  all_variants)
from ase.geometry.bravais_type_engine import identify_lattice

for lat in all_variants():
    if lat.ndim == 2:
        break

    cell = lat.tocell()

    def check(lat1):
        print('check', repr(lat), '-->', repr(lat1))
        err = np.abs(cell.cellpar() - lat1.cellpar()).max()
        assert err < 1e-5, err

    check(recognize_canonical_cell(cell))

    if lat.name not in ['MCL', 'MCLC', 'TRI']:
        stdcell, op = identify_lattice(cell, 1e-4)
        check(stdcell)
        rcell, op = cell.niggli_reduce()
        stdcell, op = identify_lattice(rcell, 1e-4)
        check(stdcell)
