import numpy as np

from ase.data import atomic_numbers
from ase.ga.utilities import closest_distances_generator
from ase.ga.bulk_utilities import CellBounds
from ase.ga.bulk_startgenerator import StartGenerator
from ase.ga.data import PrepareDB

N = 20  # number of randomly generated initial structures
blocks = [('Ag', 24)]  # the building blocks
volume = 10. * 24  # cell volume in angstrom^3
splits = {(4,): 1, (2,): 1}  # cell splitting scheme

stoichiometry = [atomic_numbers[atom] for atom, count in blocks
                 for _ in range(count)]
atom_numbers = list(set(stoichiometry))

blmin = closest_distances_generator(atom_numbers=atom_numbers,
                                    ratio_of_covalent_radii=0.5)

cellbounds = CellBounds(bounds={'phi': [0.2 * np.pi, 0.8 * np.pi],
                                'chi': [0.2 * np.pi, 0.8 * np.pi],
                                'psi': [0.2 * np.pi, 0.8 * np.pi],
                                'a': [3, 50], 'b': [3, 50], 'c': [3, 50]})

# create the starting population
sg = StartGenerator(blocks, blmin, volume, cellbounds=cellbounds,
                    splits=splits)

# create the database to store information in
da = PrepareDB(db_file_name='gadb.db',
               stoichiometry=stoichiometry)

for i in range(N):
    a = sg.get_new_candidate()
    da.add_unrelaxed_candidate(a)
