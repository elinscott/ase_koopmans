import random

from ase_koopmans import Atoms
from ase_koopmans.ga.data import PrepareDB

metals = ['Al', 'Au', 'Cu', 'Ag', 'Pd', 'Pt', 'Ni']

population_size = 10

# Create database
db = PrepareDB('fcc_alloys.db',
               population_size=population_size,
               metals=metals)

# Create starting population
for i in range(population_size):
    atoms_string = [random.choice(metals) for _ in range(4)]
    db.add_unrelaxed_candidate(Atoms(atoms_string),
                               atoms_string=''.join(atoms_string))
