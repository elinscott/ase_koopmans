import numpy as np
from ase import Atoms
from ase.ga.utilities import closest_distances_generator, atoms_too_close
from ase.ga.bulk_utilities import CellBounds
from ase.ga.bulk_startgenerator import StartGenerator
from ase.ga.bulk_crossovers import CutAndSplicePairing
from ase.ga.bulk_mutations import (SoftMutation, RotationalMutation, 
                                   StrainMutation, RattleMutation, 
                                   RattleRotationalMutation)

h2 = Atoms('H2', positions=[[0, 0, 0], [0, 0, 0.75]])
blocks = [('H', 4), ('H2O', 3), (h2, 2)]  # the building blocks
volume = 40. * sum([x[1] for x in blocks])  # cell volume in angstrom^3
splits = {(2,):1, (1,):1}  # cell splitting scheme

stoichiometry = []
for block, count in blocks:
    if type(block) == str:
        stoichiometry += list(Atoms(block).numbers) * count
    else:
        stoichiometry += list(block.numbers) * count

atom_numbers = list(set(stoichiometry))
blmin = closest_distances_generator(atom_numbers=atom_numbers,
                                    ratio_of_covalent_radii=1.3)

cellbounds = CellBounds(bounds={'phi':[0.2 * np.pi, 0.8 * np.pi],
                                'chi':[0.2 * np.pi, 0.8 * np.pi],
                                'psi':[0.2 * np.pi, 0.8 * np.pi],
                                'a':[3, 50], 'b':[3, 50], 'c':[3, 50]})

sg = StartGenerator(blocks, blmin, volume, cellbounds=cellbounds,
                    splits=splits)

# Generate 2 candidates
a1 = sg.get_new_candidate()
a1.info['confid'] = 1
a2 = sg.get_new_candidate()
a2.info['confid'] = 2

# Define and test genetic operators
pairing = CutAndSplicePairing(blmin, p1=1., p2=0., minfrac=0.15,
                              cellbounds=cellbounds, use_tags=True)

a3, desc = pairing.get_new_individual([a1, a2])
cell = a3.get_cell()
if not cellbounds.is_within_bounds(cell):
    write('to_check.traj', [a1, a2, a3])
assert cellbounds.is_within_bounds(cell)
assert not atoms_too_close(a3, blmin, use_tags=True)

strainmut = StrainMutation(blmin, stddev=0.7, cellbounds=cellbounds,
                           use_tags=True)
softmut = SoftMutation(blmin, bounds=[2., 5.], used_modes_file=None,
                       use_tags=True) 
rotmut = RotationalMutation(blmin, fraction=0.3, min_angle=0.5 * np.pi)
rattlemut = RattleMutation(blmin, rattle_prop=0.3, rattle_strength=0.5, 
                           use_tags=True)
rattlerotmut = RattleRotationalMutation(rattlemut, rotmut)

for i, mut in enumerate([strainmut, softmut, rotmut, rattlemut, rattlerotmut]):
    a = [a1, a2][i % 2]
    a3 = None
    while a3 is None:
        a3, desc = mut.get_new_individual([a])

    cell = a3.get_cell()
    assert cellbounds.is_within_bounds(cell)
    assert np.all(a3.numbers == a.numbers)
    assert not atoms_too_close(a3, blmin, use_tags=True)
