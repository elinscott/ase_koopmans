import numpy as np

from ase_koopmans.calculators.emt import EMT
from ase_koopmans.db import connect
from ase_koopmans.eos import EquationOfState
from ase_koopmans.lattice.cubic import FaceCenteredCubic

db = connect('refs.db')

metals = ['Al', 'Au', 'Cu', 'Ag', 'Pd', 'Pt', 'Ni']
for m in metals:
    atoms = FaceCenteredCubic(m)
    atoms.calc = EMT()
    e0 = atoms.get_potential_energy()
    a = atoms.cell[0][0]

    eps = 0.05
    volumes = (a * np.linspace(1 - eps, 1 + eps, 9))**3
    energies = []
    for v in volumes:
        atoms.set_cell([v**(1. / 3)] * 3, scale_atoms=True)
        energies.append(atoms.get_potential_energy())

    eos = EquationOfState(volumes, energies)
    v1, e1, B = eos.fit()

    atoms.set_cell([v1**(1. / 3)] * 3, scale_atoms=True)
    ef = atoms.get_potential_energy()

    db.write(atoms, metal=m,
             latticeconstant=v1**(1. / 3),
             energy_per_atom=ef / len(atoms))
