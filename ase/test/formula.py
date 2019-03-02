import numpy as np
from ase import Atoms

assert Atoms('MoS2').get_chemical_formula() == 'MoS2'
assert Atoms('SnO2').get_chemical_formula(mode='metal') == 'SnO2'


for sym in ['', 'Pu', 'Pu2', 'U2Pu2', 'U2(2(Pu2)H)']:
    for mode in ['all', 'reduce', 'hill', 'metal']:
        atoms = Atoms(sym)
        formula = atoms.get_chemical_formula(mode=mode)
        atoms2 = Atoms(formula)
        print(repr(sym), '->', repr(formula))
        n1 = np.sort(atoms.numbers)
        n2 = np.sort(atoms2.numbers)
        assert (n1 == n2).all()
