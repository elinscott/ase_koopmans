from random import random

from ase_koopmans.ga.cutandsplicepairing import CutAndSplicePairing
from ase_koopmans.ga.data import DataConnection
from ase_koopmans.ga.offspring_creator import OperationSelector
from ase_koopmans.ga.pbs_queue_run import PBSQueueRun
from ase_koopmans.ga.population import Population
from ase_koopmans.ga.standard_comparators import InteratomicDistanceComparator
from ase_koopmans.ga.standardmutations import (MirrorMutation,
                                               PermutationMutation,
                                               RattleMutation)
from ase_koopmans.ga.utilities import (closest_distances_generator,
                                       get_all_atom_types)
from ase_koopmans.io import write


def jtg(job_name, traj_file):
    s = '#!/bin/sh\n'
    s += '#PBS -l nodes=1:ppn=12\n'
    s += '#PBS -l walltime=48:00:00\n'
    s += '#PBS -N {0}\n'.format(job_name)
    s += '#PBS -q q12\n'
    s += 'cd $PBS_O_WORKDIR\n'
    s += 'python calc.py {0}\n'.format(traj_file)
    return s


population_size = 20
mutation_probability = 0.3

# Initialize the different components of the GA
da = DataConnection('gadb.db')
tmp_folder = 'tmp_folder/'
# The PBS queing interface is created
pbs_run = PBSQueueRun(da,
                      tmp_folder=tmp_folder,
                      job_prefix='Ag2Au2_opt',
                      n_simul=5,
                      job_template_generator=jtg)

atom_numbers_to_optimize = da.get_atom_numbers_to_optimize()
n_to_optimize = len(atom_numbers_to_optimize)
slab = da.get_slab()
all_atom_types = get_all_atom_types(slab, atom_numbers_to_optimize)
blmin = closest_distances_generator(all_atom_types,
                                    ratio_of_covalent_radii=0.7)

comp = InteratomicDistanceComparator(n_top=n_to_optimize,
                                     pair_cor_cum_diff=0.015,
                                     pair_cor_max=0.7,
                                     dE=0.02,
                                     mic=False)
pairing = CutAndSplicePairing(slab, n_to_optimize, blmin)
mutations = OperationSelector([1., 1., 1.],
                              [MirrorMutation(blmin, n_to_optimize),
                               RattleMutation(blmin, n_to_optimize),
                               PermutationMutation(n_to_optimize)])

# Relax all unrelaxed structures (e.g. the starting population)
while (da.get_number_of_unrelaxed_candidates() > 0 and
       not pbs_run.enough_jobs_running()):
    a = da.get_an_unrelaxed_candidate()
    pbs_run.relax(a)

# create the population
population = Population(data_connection=da,
                        population_size=population_size,
                        comparator=comp)

# Submit new candidates until enough are running
while (not pbs_run.enough_jobs_running() and
       len(population.get_current_population()) > 2):
    a1, a2 = population.get_two_candidates()
    a3, desc = pairing.get_new_individual([a1, a2])
    if a3 is None:
        continue
    da.add_unrelaxed_candidate(a3, description=desc)

    if random() < mutation_probability:
        a3_mut, desc = mutations.get_new_individual([a3])
        if a3_mut is not None:
            da.add_unrelaxed_step(a3_mut, desc)
            a3 = a3_mut
    pbs_run.relax(a3)

write('all_candidates.traj', da.get_all_relaxed_candidates())
