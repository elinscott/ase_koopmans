from ase.data import atomic_numbers
from ase.ga.utilities import closest_distances_generator
from ase.ga.bulk_utilities import CellBounds
from ase.ga.bulk_startgenerator import StartGenerator
from ase.ga.data import PrepareDB

# Number of random initial structures to generate
N = 20

# Target cell volume for the initial structures, in angstrom^3
volume = 240.

# Specify the 'building blocks' from which the initial structures
# will be constructed. Here we take single Ag atoms as building
# blocks, 24 in total.
blocks = [('Ag', 24)]

# From these blocks we can determine the stoichiometry
# (as a list of atomic numbers)
stoichiometry = [atomic_numbers[atom] for atom, count in blocks
                 for _ in range(count)]

# Generate a dictionary with the closest allowed interatomic distances
atom_numbers = list(set(stoichiometry))
blmin = closest_distances_generator(atom_numbers=atom_numbers,
                                    ratio_of_covalent_radii=0.5)

# Specify reasonable bounds on the minimal and maximal
# cell vector lengths (in angstrom) and angles (in degrees)
cellbounds = CellBounds(bounds={'phi': [35, 145],
                                'chi': [35, 145],
                                'psi': [35, 145],
                                'a': [3, 50], 'b': [3, 50], 'c': [3, 50]})

# Choose an (optional) 'cell splitting' scheme which basically
# controls the level of translational symmetry (within the unit
# cell) of the randomly generated structures. Here a 1:1 ratio
# of splitting factors 2 and 1 is used:
splits = {(2,): 1, (1,): 1}
# There will hence be a 50% probability that a candidate
# is constructed by repeating an randomly generated Ag12
# structure along a randomly chosen axis. In the other 50%
# of cases, no cell cell splitting will be applied.

# Initialize the random structure generator
sg = StartGenerator(blocks, blmin, volume, cellbounds=cellbounds,
                    splits=splits)

# Create the database
da = PrepareDB(db_file_name='gadb.db',
               stoichiometry=stoichiometry)

# Generate N random structures
# and add them to the database
for i in range(N):
    a = sg.get_new_candidate()
    da.add_unrelaxed_candidate(a)
