from ase.io import write
from ase.ga import get_raw_score
from ase.ga.data import DataConnection
from ase.ga.population import Population
from ase.ga.utilities import closest_distances_generator
from ase.ga.bulk_utilities import CellBounds
from ase.ga.ofp_comparator import OFPComparator
from ase.ga.offspring_creator import OperationSelector
from ase.ga.bulk_mutations import StrainMutation, SoftMutation
from ase.ga.bulk_crossovers import CutAndSplicePairing
from ga_bulk_relax import relax_one

# Connect to the database and retrieve some information
da = DataConnection('gadb.db')
atom_numbers_to_optimize = da.get_atom_numbers_to_optimize()
n_to_optimize = len(atom_numbers_to_optimize)

# Use Oganov's fingerprint functions
comp = OFPComparator(n_top=n_to_optimize, dE=1.0,
                     cos_dist_max=1e-3, rcut=10., binwidth=0.05,
                     pbc=[True, False, False], sigma=0.05, nsigma=4,
                     recalculate=False)

# Defining the cell and interatomic distance bounds that the
# candidates must obey
blmin = closest_distances_generator(atom_numbers_to_optimize, 0.5)

cellbounds = CellBounds(bounds={'phi': [35, 145],
                                'chi': [35, 145],
                                'psi': [35, 145]})

# Defining a more loose bounds for the relaxation
cb_check = CellBounds(bounds={'phi': [20, 160],
                              'chi': [20, 160],
                              'psi': [20, 160],
                              'a': [2, 60], 'b': [2, 60], 'c': [2, 60]})

# defining genetic operators:
mutation_probability = 0.6
pairing = CutAndSplicePairing(blmin, p1=1., p2=0., minfrac=0.15,
                              cellbounds=cellbounds, use_tags=False)
strainmut = StrainMutation(blmin, stddev=0.7, cellbounds=cellbounds,
                           use_tags=False)
blmin_soft = closest_distances_generator(atom_numbers_to_optimize, 0.1)

softmut = SoftMutation(blmin_soft, bounds=[2., 5.], use_tags=False)
operators = OperationSelector([4., 3., 3.],
                              [pairing, softmut, strainmut])

# relaxing the initial candidates:
while da.get_number_of_unrelaxed_candidates() > 0:
    a = da.get_an_unrelaxed_candidate()
    scaled_pos = a.get_scaled_positions(wrap=True)
    a.set_scaled_positions(scaled_pos)

    relax_one(a, cellbounds=cb_check)
    da.add_relaxed_step(a)
    if not cb_check.is_within_bounds(a.get_cell()):
        da.kill_candidate(a.info['confid'])

# create the population
population_size = 25
population = Population(data_connection=da,
                        population_size=population_size,
                        comparator=comp,
                        logfile='log.txt',
                        use_extinct=True)

# Update the scaling volume used in some operators based on a number
# of the best candidates
current_pop = population.get_current_population()
strainmut.update_scaling_volume(current_pop, w_adapt=0.5, n_adapt=4)
pairing.update_scaling_volume(current_pop, w_adapt=0.5, n_adapt=4)

# test n_to_test new candidates
n_to_test = 50
for step in range(n_to_test):
    print('Now starting configuration number {0}'.format(step))

    # Create a new candidate
    a3 = None
    while a3 is None:
        a1, a2 = population.get_two_candidates()
        a3, desc = operators.get_new_individual([a1, a2])

    # Save the unrelaxed candidate
    scaled_pos = a3.get_scaled_positions(wrap=True)
    a3.set_scaled_positions(scaled_pos)
    da.add_unrelaxed_candidate(a3, description=desc)

    # Relax the new candidate and save it
    relax_one(a3, cellbounds=cb_check)
    da.add_relaxed_step(a3)

    # If the relaxation has changed the cell parameters beyond the
    # bounds we disregard it in the population
    if not cb_check.is_within_bounds(a3.get_cell()):
        da.kill_candidate(a3.info['confid'])

    # Various updates:
    population.update()
    current_pop = population.get_current_population()

    if step % 10 == 0:
        strainmut.update_scaling_volume(current_pop, w_adapt=0.5,
                                        n_adapt=4)
        pairing.update_scaling_volume(current_pop, w_adapt=0.5, n_adapt=4)
        write('current_population.traj', current_pop)

print('GA finished after step %d' % step)
hiscore = get_raw_score(current_pop[0])
print('Highest raw score = %8.4f eV' % hiscore)

write('all_candidates.traj', da.get_all_relaxed_candidates())
write('current_population.traj', current_pop)
