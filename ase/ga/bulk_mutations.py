import numpy as np
from random import gauss
from scipy.spatial.distance import cdist
from ase import Atoms
from ase.data import chemical_symbols, covalent_radii
from ase.neighborlist import NeighborList
from ase.ga.utilities import atoms_too_close_two_sets
from ase.ga.offspring_creator import OffspringCreator
from ase.ga import standardmutations
from ase.ga.bulk_utilities import atoms_too_close, gather_atoms_by_tag,\
                                  get_rotation_matrix

class RattleMutation(standardmutations.RattleMutation):
    """ Modification of standardmutations.RattleMutation
        to allow for preserving molecular identity. """
    def __init__(self, blmin, n_top=None, rattle_strength=0.8,
                 rattle_prop=0.4, use_tags=False, verbose=False):
        standardmutations.RattleMutation.__init__(self, blmin, n_top, 
               rattle_strength=rattle_strength, rattle_prop=rattle_prop, 
               verbose=False)
        self.use_tags = use_tags

    def mutate(self, atoms):
        """ Does the actual mutation. """
        N = len(atoms) if self.n_top is None else self.n_top
        slab = atoms[:len(atoms)-N]
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
                select = np.where(tags==tag)
                cop = np.mean(pos[select], axis=0)
                if np.random.random() < self.rattle_prop:
                    ok = True
                    r = np.random.random(3)
                    pos[select] += st * (r - 0.5)

            if not ok:
                # Nothing got rattled
                continue

            top = Atoms(num, positions=pos, cell=cell, pbc=pbc, tags=tags)
            too_close = atoms_too_close(top, self.blmin, use_tags=self.use_tags)
            if not too_close:
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
            n_adapt = int(round(0.2*len(population)))
        v_new = np.mean([a.get_volume() for a in population[:n_adapt]])

        if not self.scaling_volume:
            self.scaling_volume = v_new
        else:
            volumes = [self.scaling_volume, v_new]
            weights = [1-w_adapt, w_adapt]
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
        cell = atoms.get_cell()
        vol = atoms.get_volume()
        if self.use_tags:
            tags = atoms.get_tags()
            gather_atoms_by_tag(atoms)
            pos = atoms.get_positions()

        count = 0
        invalid = True
        maxcount = 10000
        while invalid and count < maxcount:
            # generating the strain matrix:
            strain = np.identity(3)
            for i in range(3):
                for j in range(i+1):
                    if i==j:
                        strain[i,j] += gauss(0, self.stddev)
                    else:
                        epsilon = 0.5*gauss(0, self.stddev)
                        strain[i,j] += epsilon
                        strain[j,i] += epsilon

            # applying the strain:
            newcell = np.dot(strain, cell)

            # volume scaling:
            v =abs(np.linalg.det(newcell))
            if self.scaling_volume is None:
                newcell *= (vol/v)**(1./3)
            else:
                newcell *= (self.scaling_volume/v)**(1./3) 
            
            mutant = atoms.copy()
            if self.use_tags:
                transfo = np.linalg.solve(cell, newcell)
                for tag in list(set(tags)):
                    select = np.where(tags==tag)
                    cop = np.mean(pos[select], axis=0)
                    disp = np.dot(cop, transfo) - cop 
                    mutant.positions[select] += disp
                mutant.set_cell(newcell, scale_atoms=False)
            else:
                mutant.set_cell(newcell, scale_atoms=True)

            # checking constraints:
            too_close = atoms_too_close(mutant, self.blmin, 
                                        use_tags=self.use_tags)
            cell_not_ok = not self.cellbounds.is_within_bounds(newcell)
            invalid = too_close or cell_not_ok

            count += 1

        if count == maxcount:
            mutant = None
            
        return mutant


def inverse_square_model(pairs, r, parameters=None):
    ''' 
    Returns the value of the force constant k, calculated as: 
    k = k0*(r0/r)**2
    Here, the parameters k0 and r0 are taken equal to the 
    force constant and equilibrium bond length of the gas-phase
    diatomic molecule. The tabulated data have been obtained
    with PBE/def2-sv(p) in NWChem.
 
    Arguments:
    pair: Nx2 array of atomic number pairs
    r: number or array of interatomic distances (of length N)
    parameters: None (if the parameter set below is to be used),
                or dictionary with entries of the type
                (sorted atomic number pair):(k0, r0)
    '''
    if parameters is None:
        parameters = {('C', 'C'):[58.284, 1.323],
                      ('C', 'H'):[24.699, 1.15],
                      ('C', 'Mo'):[42.648, 1.685],
                      ('H', 'H'):[32.121, 0.772],
                      ('Mo', 'Mo'):[48.621, 1.982],
                      ('O', 'O'):[71.925, 1.215],
                      ('O', 'Pd'):[21.494, 1.826],
                      ('O', 'Sr'):[21.128, 1.976],
                      ('O', 'Ti'):[47.382, 1.614],
                      ('Pd', 'Pd'):[8.171, 2.499],
                      ('Sr', 'Sr'):[0.393, 4.571],
                      ('Sr', 'Ti'):[1.944, 3.22],
                      ('Ti', 'Ti'):[93.731, 1.91],
                     }

    k = []
    for i,pair in enumerate(pairs):
        key = tuple(sorted([chemical_symbols[int(j)] for j in pair]))
        k0,r0 = parameters[key]
        k.append(k0*(r0/r[i])**2) 
    return np.array(k)


class SoftMutation(OffspringCreator):
    '''
    Mutates the structure by displacing it along the lowest (nonzero)
    frequency modes found by vibrational analysis, as in:
    Lyakhov, Oganov, Valle, Comp. Phys. Comm. 181 (2010) 1623-32
    As in the reference above, the next-lowest mode is used if the
    structure has already been softmutated.
    '''
    def __init__(self, blmin, bounds=[0.5, 2.0], calculator=None,
                 fconstfunc=inverse_square_model, rcut=10.0, 
                 use_tags=False, verbose=False):
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
                    analysis. The default (None) is to use a simple
                    pairwise harmonic potential.
        fconstfunc: function that accepts a list of two atomic numbers
                    and a number of array of interatomic distances,
                    and returns the value of the force constant
                    to be used with the pairwise harmonic potential.
        rcut: cutoff radius for the pairwise harmonic potential.
        use_tags: whether to use the atomic tags to preserve
                  molecular identity.
        '''
        OffspringCreator.__init__(self, verbose)
        self.blmin = blmin
        self.bounds = bounds
        self.calc = calculator
        self.fconstfunc = fconstfunc
        self.rcut = rcut
        self.use_tags = use_tags
        self.used_modes = {}  # for storing the used modes
        self.descriptor = 'SoftMutation'
   
    def _get_pwh_hessian_(self, atoms):
        ''' Returns the Hessian matrix d2E/dxi/dxj for the pairwise
        harmonic potential. '''
        if self.use_tags:
            tags = atoms.get_tags()
        cell = atoms.get_cell()
        pos = atoms.get_positions()
        num = atoms.get_atomic_numbers()
        nat = len(atoms) 

        # build neighborlist
        nl = NeighborList([self.rcut/2.]*nat, skin=0., bothways=True,
                          self_interaction=False)
        nl.build(atoms)

        # constructing the hessian
        hessian = np.zeros((nat*3, nat*3))
        large_fc = 1e5  # high force constant for same-tag atoms
        for i in range(nat):
            indices, offsets = nl.get_neighbors(i)
            for j in range(nat):
                if i==j:
                    p = pos[indices] + np.dot(offsets, cell)
                    r = cdist(p, [pos[i]])
                    v = p - pos[i]
                    pairs = np.vstack((num[indices], np.zeros(len(indices)))).T 
                    pairs[:,1] = num[i]
                    fc = self.fconstfunc(pairs, r[:,0])

                    if self.use_tags:
                        tag_indices = np.where(tags==tags[i])[0]
                        for tag_index in tag_indices:
                            if tag_index==i:
                                continue
                            select = np.where(indices==tag_index)
                            fc[select] = large_fc 
                                    
                    v /= r
                    for k in range(3):
                        for l in range(3):
                            index1 = 3*i+k
                            index2 = 3*j+l
                            h = np.sum(np.dot(fc*v[:,k], v[:,l]))
                            hessian[index1, index2] = h
                else:
                    for m,index in enumerate(indices):
                        if index != j: 
                            continue
                        v = pos[index] + np.dot(offsets[m], cell) - pos[i]
                        r = np.linalg.norm(v)
                        v /= r
                        for k in range(3):
                            for l in range(3):
                                index1 = 3*i+k
                                index2 = 3*j+l
                                fc = self.fconstfunc([[num[i], num[j]]], [r])
                                if self.use_tags and tags[i]==tags[j]:
                                    fc = large_fc
                                h = -1*fc*v[k]*v[l]
                                hessian[index1, index2] += h
        return hessian

    def _get_hessian_(self, atoms, dx):
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
        hessian = np.zeros((3*nat, 3*nat))
        for i in range(3*nat):
            row = np.zeros(3*nat)
            for direction in [-1, 1]:
                disp = np.zeros(3)
                disp[i%3] = direction*dx
                pos_disp = np.copy(pos)
                pos_disp[i/3] += disp
                atoms.positions = pos_disp
                f = atoms.get_forces()
                row += -1*direction*f.flatten()
            row /= (2.*dx)
            hessian[i] = row
        hessian += np.copy(hessian).T
        hessian *= 0.5
        return hessian

    def _calculate_normal_modes_(self, atoms, dx=0.02, massweighing=False):
        '''Performs the vibrational analysis.'''
        if self.calc is None:
            hessian = self._get_pwh_hessian_(atoms)
        else:
            hessian = self._get_hessian_(atoms,dx)

        if massweighing:
            m = np.array([np.repeat(atoms.get_masses()**-0.5,3)])
            hessian *= (m*m.T)

        eigvals,eigvecs = np.linalg.eigh(hessian)
        modes = {eigval:eigvecs[:,i] for i,eigval in enumerate(eigvals)}
        return modes

    def _animate_mode_(self, atoms, mode, nim=30, amplitude=1.0):
        '''Returns an Atoms object showing an animation of the mode.'''
        pos = atoms.get_positions()
        mode = mode.reshape(np.shape(pos))
        animation = []
        for i in range(nim):
            newpos = pos + amplitude*mode*np.sin(i*2*np.pi/nim)
            image = atoms.copy()
            image.positions = newpos
            animation.append(image)	
        return animation

    def get_new_individual(self, parents):
        f = parents[0]

        indi = self.mutate(f)
        if indi is None:
            return indi, 'mutation: soft'
            
        indi = self.initialize_individual(f, indi)
        indi.info['data']['parents'] = [f.info['confid']]

        return self.finalize_individual(indi), 'mutation: soft'

    def mutate(self, a):
        """ Does the actual mutation. """
        pos = atoms.get_positions()
        num = atoms.get_atomic_numbers()        

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
        while amplitude*largest_norm < self.bounds[1]:
            newpos = pos + direction*amplitude*mode
            mutant.set_positions(newpos)
            mutant.wrap()
            too_close = atoms_too_close(mutant, self.blmin, 
                                        use_tags=self.use_tags)
            if too_close:
                amplitude -= increment
                newpos = pos + direction*amplitude*mode
                mutant.set_positions(newpos)
                mutant.wrap()
                break

            if direction == 1:
                direction = -1
            else:
                direction = 1
                amplitude += increment

        if amplitude*largest_norm < self.bounds[0]:
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
    arXiv:1204.4756v2
    """
    def __init__(self, blmin, n_top=None, probability=0.33, tags=None, 
                 test_dist_to_slab=True, verbose=False):
        """ Parameters:
        blmin: closest allowed distances
        n_top: number of atoms to optimize; if None, all are included.
        probability: probability with which a moiety is rotated.
        tags: None or list of integers, specify respectively whether 
              all moieties or only those with matching tags are 
              eligible for rotation.
        test_dist_to_slab: whether also the distances to the slab
                           should be checked to satisfy the blmin.
        """
        OffspringCreator.__init__(self, verbose)
        self.blmin = blmin
        self.n_top = n_top
        self.probability = probability
        self.tags = tags
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
        slab = atoms[:len(atoms)-N]
        atoms = atoms[-N:]

        mutant = atoms.copy()
        gather_atoms_by_tag(mutant)
        pos = mutant.get_positions() 
        tags = mutant.get_tags()
        eligible_tags = tags if self.tags is None else self.tags

        indices = {}
        for tag in list(set(tags)):
            hits = np.where(tags==tag)[0]
            if len(hits) > 1 and tag in eligible_tags:
                indices[tag] = hits

        n_rot = int(np.ceil(len(indices)*self.probability))
        chosen_tags = np.random.choice(indices.keys(), size=n_rot,
                                       replace=False)

        too_close = True
        count = 0
        maxcount = 10000
        while too_close and count < maxcount:
            newpos = np.copy(pos)
            for tag in chosen_tags:
                p = np.copy(newpos[indices[tag]])
                cop = np.mean(p, axis=0)
                axis = np.random.random(3)
                axis /= np.linalg.norm(axis)
                angle = 2*np.pi*np.random.random()
                m = get_rotation_matrix(axis, angle)
                newpos[indices[tag]] = np.dot(m, (p-cop).T).T + cop

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
