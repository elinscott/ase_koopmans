"""Mutation operations intended for bulk structures.
If you find this implementation useful in your work,
please cite:
    M. Van den Bossche, Henrik Gronbeck, B. Hammer,
    J. Chem. Theory Comput., doi:10.1021/acs.jctc.8b00039
in addition to the papers mentioned in the docstrings."""
import json
import numpy as np
from random import gauss
from ase import Atoms
from ase.data import covalent_radii
from ase.neighborlist import NeighborList
from ase.build import niggli_reduce
from ase.ga.offspring_creator import OffspringCreator
from ase.ga import standardmutations
from ase.ga.utilities import (atoms_too_close, atoms_too_close_two_sets,
                              gather_atoms_by_tag)
from ase.ga.bulk_utilities import get_rotation_matrix
try:
    from scipy.spatial.distance import cdist
    have_scipy = True
except ImportError
    have_scipy = False


class RattleMutation(standardmutations.RattleMutation):
    """ Modification of standardmutations.RattleMutation
        to allow for preserving molecular identity. """

    def __init__(self, blmin, n_top=None, rattle_strength=0.8,
                 rattle_prop=0.4, use_tags=False, test_dist_to_slab=True,
                 verbose=False):
        standardmutations.RattleMutation.__init__(self, blmin, n_top,
                                                  rattle_strength=rattle_strength, 
                                                  rattle_prop=rattle_prop,
                                                  verbose=verbose)
        self.use_tags = use_tags
        self.test_dist_to_slab = test_dist_to_slab

    def mutate(self, atoms):
        """ Does the actual mutation. """
        N = len(atoms) if self.n_top is None else self.n_top
        slab = atoms[:len(atoms) - N]
        atoms = atoms[-N:]
        tags = atoms.get_tags() if self.use_tags else np.arange(N)
        pos_ref = atoms.get_positions()
        num = atoms.get_atomic_numbers()
        cell = atoms.get_cell()
        pbc = atoms.get_pbc()
        st = 2. * self.rattle_strength

        count = 0
        maxcount = 1000
        too_close = True
        while too_close and count < maxcount:
            count += 1
            pos = pos_ref.copy()
            ok = False
            for tag in list(set(tags)):
                select = np.where(tags == tag)
                if np.random.random() < self.rattle_prop:
                    ok = True
                    r = np.random.random(3)
                    pos[select] += st * (r - 0.5)

            if not ok:
                # Nothing got rattled
                continue

            top = Atoms(num, positions=pos, cell=cell, pbc=pbc, tags=tags)
            too_close = atoms_too_close(
                top, self.blmin, use_tags=self.use_tags)
            if not too_close and self.test_dist_to_slab:
                too_close = atoms_too_close_two_sets(top, slab, self.blmin)

        if count == maxcount:
            return None

        mutant = slab + top
        return mutant


class PermutationMutation(standardmutations.PermutationMutation):
    """ Modification of standardmutations.PermutationMutation
        to allow for preserving molecular identity. """

    def __init__(self, blmin, n_top=None, probability=0.33, use_tags=False,
                 test_dist_to_slab=True, verbose=False):
        standardmutations.PermutationMutation.__init__(self, n_top,
                                                       probability=probability, verbose=verbose)
        self.blmin = blmin
        self.use_tags = use_tags
        self.test_dist_to_slab = test_dist_to_slab

    def mutate(self, atoms):
        """ Does the actual mutation. """
        N = len(atoms) if self.n_top is None else self.n_top
        slab = atoms[:len(atoms) - N]
        atoms = atoms[-N:]
        if self.use_tags:
            gather_atoms_by_tag(atoms)
        tags = atoms.get_tags() if self.use_tags else np.arange(N)
        pos_ref = atoms.get_positions()
        num = atoms.get_atomic_numbers()
        cell = atoms.get_cell()
        pbc = atoms.get_pbc()
        symbols = atoms.get_chemical_symbols()

        unique_tags = list(set(tags))
        n = len(unique_tags)
        swaps = int(np.ceil(n * self.probability / 2.))

        sym = []
        for tag in list(set(unique_tags)):
            indices = np.where(tags == tag)[0]
            s = ''.join([symbols[j] for j in indices])
            sym.append(s)
        assert len(list(set(sym))) > 1

        count = 0
        maxcount = 1000
        too_close = True
        while too_close and count < maxcount:
            count += 1
            pos = pos_ref.copy()
            for _ in range(swaps):
                i = j = 0
                while sym[i] == sym[j]:
                    i = np.random.randint(0, high=n)
                    j = np.random.randint(0, high=n)
                ind1 = np.where(tags == i)
                ind2 = np.where(tags == j)
                cop1 = np.mean(pos[ind1], axis=0)
                cop2 = np.mean(pos[ind2], axis=0)
                pos[ind1] += cop2 - cop1
                pos[ind2] += cop1 - cop2

            top = Atoms(num, positions=pos, cell=cell, pbc=pbc, tags=tags)
            too_close = atoms_too_close(
                top, self.blmin, use_tags=self.use_tags)
            if not too_close and self.test_dist_to_slab:
                too_close = atoms_too_close_two_sets(top, slab, self.blmin)

        if count == maxcount:
            return None

        mutant = slab + top
        return mutant


class PermuStrainMutation(OffspringCreator):
    """ Combination of PermutationMutation and StrainMutation, see also:
    Lonie, Zurek, Comp. Phys. Comm. 182 (2011) 372-387
    """

    def __init__(self, permutationmutation, strainmutation, verbose=False):
        """
        permutationmutation: instance of a mutation that permutes 
                             atom types 
        strainmutation: instance of a mutation that mutates by straining
        """
        OffspringCreator.__init__(self, verbose)
        self.permutationmutation = permutationmutation
        self.strainmutation = strainmutation

    def get_new_individual(self, parents):
        f = parents[0]

        indi = self.mutate(f)
        if indi is None:
            return indi, 'mutation: permustrain'

        indi = self.initialize_individual(f, indi)
        indi.info['data']['parents'] = [f.info['confid']]

        return self.finalize_individual(indi), 'mutation: permustrain'

    def mutate(self, atoms):
        """ Does the actual mutation. """
        mutant = self.permutationmutation.mutate(atoms)
        if mutant is not None:
            mutant = self.strainmutation.mutate(mutant)
        return mutant


class StrainMutation(OffspringCreator):
    """ Mutates a candidate by applying a randomly generated strain.
    See also:
    Lonie, Zurek, Comp. Phys. Comm. 182 (2011) 372-387
    Glass, Oganov, Hansen, Comp. Phys. Comm. 175 (2006) 713-720

    After initialization of the mutation, a scaling volume
    (to which each mutated structure is scaled before checking the
    constraints) is typically generated from the population, 
    which is then also occasionally updated in the course of the
    GA run.
    """

    def __init__(self, blmin, cellbounds=None, stddev=0.7, use_tags=False,
                 verbose=False):
        """ Parameters:
        blmin: dict with the minimal interatomic distances
        cellbounds: ase.ga.bulk_utilities.CellBounds instance 
                    describing limits on the cell shape
        stddev: standard deviation used in the generation of the strain
                matrix elements
        use_tags: whether to use the atomic tags to preserve
                  molecular identity.
        """
        OffspringCreator.__init__(self, verbose)
        self.blmin = blmin
        self.cellbounds = cellbounds
        self.stddev = stddev
        self.use_tags = use_tags
        self.scaling_volume = None
        self.descriptor = 'StrainMutation'
        self.min_inputs = 1

    def update_scaling_volume(self, population, w_adapt=0.5, n_adapt=0):
        """Function to initialize or update the scaling volume in a GA run."""
        if not n_adapt:
            # if not set, take best 20% of the population
            n_adapt = int(round(0.2 * len(population)))
        v_new = np.mean([a.get_volume() for a in population[:n_adapt]])

        if not self.scaling_volume:
            self.scaling_volume = v_new
        else:
            volumes = [self.scaling_volume, v_new]
            weights = [1 - w_adapt, w_adapt]
            self.scaling_volume = np.average(volumes, weights=weights)

    def get_new_individual(self, parents):
        f = parents[0]

        indi = self.mutate(f)
        if indi is None:
            return indi, 'mutation: strain'

        indi = self.initialize_individual(f, indi)
        indi.info['data']['parents'] = [f.info['confid']]

        return self.finalize_individual(indi), 'mutation: strain'

    def mutate(self, atoms):
        """ Does the actual mutation. """
        cell_ref = atoms.get_cell()
        pos_ref = atoms.get_positions()
        vol = atoms.get_volume()
        if self.use_tags:
            tags = atoms.get_tags()
            gather_atoms_by_tag(atoms)
            pos = atoms.get_positions()

        mutant = atoms.copy()
        if self.cellbounds is not None:
            if not self.cellbounds.is_within_bounds(cell_ref):
                niggli_reduce(mutant)

        count = 0
        too_close = True
        maxcount = 1000
        while too_close and count < maxcount:
            mutant.set_cell(cell_ref, scale_atoms=False)
            mutant.set_positions(pos_ref)

            # generating the strain matrix:
            strain = np.identity(3)
            for i in range(3):
                for j in range(i + 1):
                    if i == j:
                        strain[i, j] += gauss(0, self.stddev)
                    else:
                        epsilon = 0.5 * gauss(0, self.stddev)
                        strain[i, j] += epsilon
                        strain[j, i] += epsilon

            # applying the strain:
            cell_new = np.dot(strain, cell_ref)

            # volume scaling:
            v = abs(np.linalg.det(cell_new))
            if self.scaling_volume is None:
                cell_new *= (vol / v)**(1. / 3)
            else:
                cell_new *= (self.scaling_volume / v)**(1. / 3)

            # check cell dimensions:
            if not self.cellbounds.is_within_bounds(cell_new):
                continue

            if self.use_tags:
                transfo = np.linalg.solve(cell_ref, cell_new)
                for tag in list(set(tags)):
                    select = np.where(tags == tag)
                    cop = np.mean(pos[select], axis=0)
                    disp = np.dot(cop, transfo) - cop
                    mutant.positions[select] += disp

            mutant.set_cell(cell_new, scale_atoms=not self.use_tags)

            # check distances:
            too_close = atoms_too_close(mutant, self.blmin,
                                        use_tags=self.use_tags)
            count += 1

        if count == maxcount:
            mutant = None

        return mutant


def get_number_of_valence_electrons(Z):
    ''' Return the number of valence electrons for the element with
        atomic number Z, simply based on its periodic table group '''
    groups = [[], [1, 3, 11, 19, 37, 55, 87], [2, 4, 12, 20, 38, 56, 88],
              [21, 39, 57, 89]]

    for i in range(9):
        groups.append(i + np.array([22, 40, 72, 104]))

    for i in range(6):
        groups.append(i + np.array([5, 13, 31, 49, 81, 113]))

    for i, group in enumerate(groups):
        if Z in group:
            nval = i if i < 13 else i - 10
            break
    else:
        raise ValueError('Z=%d not included in this dataset.' % Z) 

    return nval


class SoftMutation(OffspringCreator):
    '''
    Mutates the structure by displacing it along the lowest (nonzero)
    frequency modes found by vibrational analysis, as in:
    Lyakhov, Oganov, Valle, Comp. Phys. Comm. 181 (2010) 1623-32
    As in the reference above, the next-lowest mode is used if the
    structure has already been softmutated.
    '''

    def __init__(self, blmin, bounds=[0.5, 2.0], calculator=None, rcut=10.,
                 used_modes_file='used_modes.json', use_tags=False,
                 verbose=False):
        '''
        blmin: dictionary with closest allowed interatomic distances.
        bounds: lower and upper limits (in Angstrom) for the largest 
                atomic displacement in the structure. For a given mode,  
                the algorithm starts at zero amplitude and increases 
                it until either blmin is violated or the largest 
                displacement exceeds the provided upper bound).
                If the largest displacement in the resulting structure
                is lower than the provided lower bound, the mutant is
                considered too similar to the parent and None is 
                returned.
        calculator: the calculator to be used in the vibrational 
                    analysis. The default (None) is to determine the
                    the force constants via the "bond electronegativity"
                    model described in the reference above.
        rcut: cutoff radius for the pairwise harmonic potential.
        used_modes_file: name of json dump file where previously used 
                  modes will be stored (and read). If None, no such 
                  file will be used.
        use_tags: whether to use the atomic tags to preserve
                  molecular identity.
        '''
        OffspringCreator.__init__(self, verbose)
        self.blmin = blmin
        self.bounds = bounds
        self.calc = calculator
        self.rcut = rcut
        self.used_modes_file = used_modes_file
        self.use_tags = use_tags
        self.descriptor = 'SoftMutation'

        assert have_scipy, 'SoftMutation requires SciPy'

        self.used_modes = {}
        if self.used_modes_file is not None:
            try:
                self.read_used_modes(self.used_modes_file)
            except IOError:
                # file doesn't exist (yet)
                pass

    def _get_bem_hessian_(self, atoms):
        ''' Returns the Hessian matrix d2E/dxi/dxj using the bond 
            electronegativity model '''
        if self.use_tags:
            tags = atoms.get_tags()
        cell = atoms.get_cell()
        pos = atoms.get_positions()
        num = atoms.get_atomic_numbers()
        nat = len(atoms)

        # build neighborlist
        nl = NeighborList([self.rcut / 2.] * nat, skin=0., bothways=True,
                          self_interaction=False)
        nl.update(atoms)

        # computing the force constants
        s_norms = []
        valence_states = []
        r_cov = []
        for i in range(nat):
            indices, offsets = nl.get_neighbors(i)
            p = pos[indices] + np.dot(offsets, cell)
            r = cdist(p, [pos[i]])
            r_ci = covalent_radii[num[i]]
            s = 0.
            for j in indices:
                d = r[j] - r_ci - covalent_radii[num[j]]
                s += np.exp(-d / 0.37)
            s_norms.append(s)
            valence_states.append(get_number_of_valence_electrons(num[i]))
            r_cov.append(r_ci)

        fconst = [] 
        for i in range(nat):
            indices, offsets = nl.get_neighbors(i)
            p = pos[indices] + np.dot(offsets, cell)
            r = cdist(p, [pos[i]])
            n = num[indices]
            fc = []
            for j in indices:
                d = r[j] - r_cov[i] - r_cov[j]
                chi_ik = 0.481 * valence_states[num[i]] / (r_cov[i] + 0.5 * d)
                chi_jk = 0.481 * valence_states[num[j]] / (r_cov[j] + 0.5 * d)
                cn_ik = np.exp(-d / 0.37) / s_norms[i]
                cn_jk = np.exp(-d / 0.37) / s_norms[j]
                fc.append(np.sqrt(chi_ik * chi_jk / (cn_ik * cn_jk)))
            fconst.append(fc)

        # constructing the hessian
        hessian = np.zeros((nat * 3, nat * 3))
        large_fc = 1e5  # high force constant for same-tag atoms
        for i in range(nat):
            indices, offsets = nl.get_neighbors(i)
            fcs = np.array(fconst[i])
            for j in range(nat):
                if i == j:
                    p = pos[indices] + np.dot(offsets, cell)
                    r = cdist(p, [pos[i]])
                    v = p - pos[i]
                    fc = fcs.copy()
  
                    if self.use_tags:
                        tag_indices = np.where(tags == tags[i])[0]
                        for tag_index in tag_indices:
                            if tag_index == i:
                                continue
                            select = np.where(indices == tag_index)
                            fc[select] = large_fc

                    v /= r
                    for k in range(3):
                        for l in range(3):
                            index1 = 3 * i + k
                            index2 = 3 * j + l
                            h = np.sum(np.dot(fc * v[:, k], v[:, l]))
                            hessian[index1, index2] = h
                else:
                    for m, index in enumerate(indices):
                        if index != j:
                            continue
                        v = pos[index] + np.dot(offsets[m], cell) - pos[i]
                        r = np.linalg.norm(v)
                        v /= r
                        for k in range(3):
                            for l in range(3):
                                index1 = 3 * i + k
                                index2 = 3 * j + l
                                fc = fcs[m] 
                                if self.use_tags and tags[i] == tags[j]:
                                    fc = large_fc
                                h = -1 * fc * v[k] * v[l]
                                hessian[index1, index2] += h
        return hessian

    def _get_calc_hessian_(self, atoms, dx):
        ''' 
        Returns the Hessian matrix d2E/dxi/dxj using self.calc as
        calculator, through a first-order central difference scheme with
        displacements dx. In principle we could use the ase.vibrations 
        module, but that one involves alot of I/O operations, so this is
        a more light-weight version better suited for a GA mutation.
        '''
        pos = atoms.get_positions()
        atoms.set_calculator(self.calc)
        nat = len(atoms)
        hessian = np.zeros((3 * nat, 3 * nat))
        for i in range(3 * nat):
            row = np.zeros(3 * nat)
            for direction in [-1, 1]:
                disp = np.zeros(3)
                disp[i % 3] = direction * dx
                pos_disp = np.copy(pos)
                pos_disp[i / 3] += disp
                atoms.positions = pos_disp
                f = atoms.get_forces()
                row += -1 * direction * f.flatten()
            row /= (2. * dx)
            hessian[i] = row
        hessian += np.copy(hessian).T
        hessian *= 0.5
        return hessian

    def _calculate_normal_modes_(self, atoms, dx=0.02, massweighing=False):
        '''Performs the vibrational analysis.'''
        if self.calc is None:
            hessian = self._get_bem_hessian_(atoms)
        else:
            hessian = self._get_calc_hessian_(atoms, dx)

        if massweighing:
            m = np.array([np.repeat(atoms.get_masses()**-0.5, 3)])
            hessian *= (m * m.T)

        eigvals, eigvecs = np.linalg.eigh(hessian)
        modes = {eigval: eigvecs[:, i] for i, eigval in enumerate(eigvals)}
        return modes

    def _animate_mode_(self, atoms, mode, nim=30, amplitude=1.0):
        '''Returns an Atoms object showing an animation of the mode.'''
        pos = atoms.get_positions()
        mode = mode.reshape(np.shape(pos))
        animation = []
        for i in range(nim):
            newpos = pos + amplitude * mode * np.sin(i * 2 * np.pi / nim)
            image = atoms.copy()
            image.positions = newpos
            animation.append(image)
        return animation

    def read_used_modes(self, filename):
        ''' Read used modes from json file. '''
        with open(filename, 'r') as f:
            modes = json.load(f)
            self.used_modes = {int(k): modes[k] for k in modes}
        return

    def write_used_modes(self, filename):
        ''' Dump used modes to json file. '''
        with open(filename, 'w') as f:
            json.dump(self.used_modes, f)
        return

    def get_new_individual(self, parents):
        f = parents[0]

        indi = self.mutate(f)
        if indi is None:
            return indi, 'mutation: soft'

        indi = self.initialize_individual(f, indi)
        indi.info['data']['parents'] = [f.info['confid']]

        return self.finalize_individual(indi), 'mutation: soft'

    def mutate(self, atoms):
        """ Does the actual mutation. """
        pos = atoms.get_positions()
        modes = self._calculate_normal_modes_(atoms)

        # Select the mode along which we want to move the atoms;
        # The first 3 translational modes as well as previously
        # applied modes are discarded.

        keys = np.array(sorted(modes))
        index = 3
        confid = atoms.info['confid']
        if confid in self.used_modes:
            while index in self.used_modes[confid]:
                index += 1
            self.used_modes[confid].append(index)
        else:
            self.used_modes[confid] = [index]

        if self.used_modes_file is not None:
            self.write_used_modes(self.used_modes_file)

        key = keys[index]
        mode = modes[key].reshape(np.shape(pos))

        # Find a suitable amplitude for translation along the mode;
        # at every trial amplitude both positive and negative
        # directions are tried.

        mutant = atoms.copy()
        amplitude = 0.
        increment = 0.1
        direction = 1
        largest_norm = np.max(np.apply_along_axis(np.linalg.norm, 1, mode))
        while amplitude * largest_norm < self.bounds[1]:
            newpos = pos + direction * amplitude * mode
            mutant.set_positions(newpos)
            mutant.wrap()
            too_close = atoms_too_close(mutant, self.blmin,
                                        use_tags=self.use_tags)
            if too_close:
                amplitude -= increment
                newpos = pos + direction * amplitude * mode
                mutant.set_positions(newpos)
                mutant.wrap()
                break

            if direction == 1:
                direction = -1
            else:
                direction = 1
                amplitude += increment

        if amplitude * largest_norm < self.bounds[0]:
            mutant = None

        return mutant


class RotationalMutation(OffspringCreator):
    """ Mutates a candidate by applying random rotations 
    to multi-atom moieties in the structure (atoms with
    the same tag are considered part of one such moiety).
    Only performs whole-molecule rotations, no internal
    rotations.   

    See also:
    Zhu Q., Oganov A.R., Glass C.W., Stokes H.T,
      Acta Cryst. (2012), B68, 215-226.
    """

    def __init__(self, blmin, n_top=None, fraction=0.33, tags=None,
                 min_angle=1.57, test_dist_to_slab=True, verbose=False):
        """ Parameters:
        blmin: closest allowed distances
        n_top: number of atoms to optimize; if None, all are included.
        fraction: fraction of the moieties to be rotated.
        tags: None or list of integers, specify respectively whether 
              all moieties or only those with matching tags are 
              eligible for rotation.
        min_angle: minimal angle (in radians) for each rotation;
                   should lie in the interval [0, pi].
        test_dist_to_slab: whether also the distances to the slab
                           should be checked to satisfy the blmin.
        """
        OffspringCreator.__init__(self, verbose)
        self.blmin = blmin
        self.n_top = n_top
        self.fraction = fraction
        self.tags = tags
        self.min_angle = min_angle
        self.test_dist_to_slab = test_dist_to_slab
        self.descriptor = 'RotationalMutation'
        self.min_inputs = 1

    def get_new_individual(self, parents):
        f = parents[0]

        indi = self.mutate(f)
        if indi is None:
            return indi, 'mutation: rotational'

        indi = self.initialize_individual(f, indi)
        indi.info['data']['parents'] = [f.info['confid']]

        return self.finalize_individual(indi), 'mutation: rotational'

    def mutate(self, atoms):
        """ Does the actual mutation. """
        N = len(atoms) if self.n_top is None else self.n_top
        slab = atoms[:len(atoms) - N]
        atoms = atoms[-N:]

        mutant = atoms.copy()
        gather_atoms_by_tag(mutant)
        pos = mutant.get_positions()
        tags = mutant.get_tags()
        eligible_tags = tags if self.tags is None else self.tags

        indices = {}
        for tag in list(set(tags)):
            hits = np.where(tags == tag)[0]
            if len(hits) > 1 and tag in eligible_tags:
                indices[tag] = hits

        n_rot = int(np.ceil(len(indices) * self.fraction))
        chosen_tags = np.random.choice(list(indices.keys()), size=n_rot,
                                       replace=False)

        too_close = True
        count = 0
        maxcount = 10000
        while too_close and count < maxcount:
            newpos = np.copy(pos)
            for tag in chosen_tags:
                p = np.copy(newpos[indices[tag]])
                cop = np.mean(p, axis=0)

                if len(p) == 2:
                    line = (p[1] - p[0]) / np.linalg.norm(p[1] - p[0])
                    while True:
                        axis = np.random.random(3)
                        axis /= np.linalg.norm(axis)
                        a = np.arccos(np.dot(axis, line))
                        if np.pi / 4 < a < np.pi * 3 / 4:
                            break
                else:
                    axis = np.random.random(3)
                    axis /= np.linalg.norm(axis)

                angle = self.min_angle
                angle += 2 * (np.pi - self.min_angle) * np.random.random()

                m = get_rotation_matrix(axis, angle)
                newpos[indices[tag]] = np.dot(m, (p - cop).T).T + cop

            mutant.set_positions(newpos)
            mutant.wrap()
            too_close = atoms_too_close(mutant, self.blmin, use_tags=True)
            count += 1

            if not too_close and self.test_dist_to_slab:
                too_close = atoms_too_close_two_sets(slab, mutant, self.blmin)

        if count == maxcount:
            mutant = None
        else:
            mutant = slab + mutant

        return mutant


class RattleRotationalMutation(OffspringCreator):
    """ Combination of Rattle and RotationalMutations """

    def __init__(self, rattlemutation, rotationalmutation, verbose=False):
        """
        rattlemutation: instance that rattles atoms
        rotationalmutation: instance of a mutation that rotates moieties
        """
        OffspringCreator.__init__(self, verbose)
        self.rattlemutation = rattlemutation
        self.rotationalmutation = rotationalmutation
        self.descriptor = 'rattlerotational'

    def get_new_individual(self, parents):
        f = parents[0]

        indi = self.mutate(f)
        if indi is None:
            return indi, 'mutation rattlerotational'

        indi = self.initialize_individual(f, indi)
        indi.info['data']['parents'] = [f.info['confid']]

        return self.finalize_individual(indi), 'mutation: rattlerotational'

    def mutate(self, atoms):
        """ Does the actual mutation. """
        mutant = self.rattlemutation.mutate(atoms)
        if mutant is not None:
            mutant = self.rotationalmutation.mutate(mutant)
        return mutant
