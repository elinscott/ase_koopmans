import numpy as np
from random import random, randrange
from ase import Atoms
from ase.ga.bulk_utilities import atoms_too_close
from ase.ga.offspring_creator import OffspringCreator

class Position(object):
    """Position helper object.

    This class is just a simple representation used by the pairing
    operator.

    Parameters:

    position: [x, y, z] coordinate
    number: Atomic species at the position
    distance: Signed distance to the cutting plane
    origin: Either 0 or 1 and determines which side of the plane
    the position should be at.

    """
    def __init__(self, position, scaled_position, number, distance, origin):
        self.number = number
        self.position = position
        self.scaled_position = scaled_position
        self.distance = distance
        self.origin = origin

    def to_use(self):
        """ Method which tells if this position is at the right side.
        """
        if self.distance > 0. and self.origin == 0:
            return True
        elif self.distance < 0. and self.origin == 1:
            return True
        else:
            return False


class CutAndSplicePairing(OffspringCreator):
    """ Parameters:
    n_top   : The number of atoms to optimize
    blmin   : Dictionary with pairs of atom numbers and the closest
              distance these can have to each other.
    p1      : probability that a parent is shifted over a random
              distance along the normal of the cutting plane 
    p2      : same as p1, but for shifting along the two directions 
              in the cutting plane
    minfrac : minimal fraction of atoms a parent must contribute 
              to the child
    use_tags: whether to use the atomic tags to preserve
              molecular identity.
 
    For more information, see e.g.
    Glass, Oganov, Hansen, Comp. Phys. Comm. 175 (2006) 713-720
    Lonie, Zurek, Comp. Phys. Comm. 182 (2011) 372-387
    """
    def __init__(self, n_top, blmin, p1=1., p2=0.05, minfrac=None,  
                 use_tags=False, verbose=False):
        OffspringCreator.__init__(self, verbose)
        self.blmin = blmin
        self.n_top = n_top
        self.p1 = p1
        self.p2 = p2
        self.minfrac = minfrac
        self.scaling_volume = None
        self.use_tags = use_tags
        self.descriptor = 'CutAndSplicePairing'
        self.min_inputs = 2

    def update_scaling_volume(self,population, w_adapt=0.5, n_adapt=0):
        ''' Updates the scaling volume that is used in the pairing
        w_adapt: weight of the new vs the old scaling volume
        n_adapt: number of best candidates in the population that 
                 are used to calculate the new scaling volume 
        '''
        if not n_adapt:
            # take best 20% of the population
            n_adapt = int(round(0.2*len(population)))
        v_new = np.mean([ a.get_volume() for a in population[:n_adapt] ])

        if not self.scaling_volume:
            self.scaling_volume = v_new
        else:
            volumes = [self.scaling_volume, v_new]
            weights = [1-w_adapt, w_adapt]
            self.scaling_volume = np.average(volumes, weights=weights)

    def _get_pairing_(self, a1, a2, direction=None, fraction=None):
        """ 
        Creates a child from two parents using the given cutting plane
        Does not check whether atoms are too close. 
        direction: direction of the cutting surface normal (0, 1 or 2) 
        fraction: fraction of the lattice vector along which  
                  the cut is made 
        """
        N = len(a1)
        num = a1.numbers[:]
        unique_types = list(set(num))

        types = dict()
        for u in unique_types:
            types[u] = sum(num == u)

        # Generate list of all atoms
        cell1 = a1.get_cell()
        cell2 = a2.get_cell()
        scalpos1 = a1.get_scaled_positions()
        scalpos2 = a2.get_scaled_positions()
        p1 = [Position(a.position, scalpos1[i], a.number,
            scalpos1[i,direction]-fraction, origin=0) for i,a in enumerate(a1)]
        p2 = [Position(a.position, scalpos2[i], a.number, 
            scalpos2[i,direction]-fraction, origin=1) for i,a in enumerate(a2)]

        all_points = p1
        all_points.extend(p2)

        # Sort these by their atomic number
        all_points.sort(key=lambda x: x.number, reverse=True)

        # For each atom type make the pairing
        unique_types.sort()
        use_total = dict()
        for u in unique_types:
            used = []
            not_used = []
            # The list is looked trough in
            # reverse order so atoms can be removed
            # from the list along the way.
            for i in reversed(range(len(all_points))):
                # If there are no more atoms of this type
                if all_points[i].number != u:
                    break
                # Check if the atom should be included
                if all_points[i].to_use():
                    used.append(all_points.pop(i))
                else:
                    not_used.append(all_points.pop(i))

            assert len(used) + len(not_used) == types[u] * 2

            # While we have too few of the given atom type
            while len(used) < types[u]:
                r = random()
                # origin = 0 => provides atoms if pos > fraction
                # origin = 1 => provides atoms if pos < fraction
                pick = 0 if r > fraction else 1
                interval = [0,fraction] if r < fraction else [fraction,1]
                indices = []
                for index,point in enumerate(not_used):
                    cond1 = interval[0] <= point.scaled_position[direction]
                    cond2 = point.scaled_position[direction] <= interval[1]
                    if cond1 and cond2:
                        if point.origin != pick:
                            indices.append(index)
                choice = randrange(0,len(indices))
                used.append(not_used.pop(choice))

            # While we have too many of the given atom type
            while len(used) > types[u]:
                # remove randomly:
                index = randrange(0,len(used))
                not_used.append(used.pop(index))

            use_total[u] = used

        n_tot = sum([len(ll) for ll in use_total.values()])
        assert n_tot == N
 
        # pair the cells:
        r = random()
        newcell = np.average([cell1,cell2], weights=[r,1-r], axis=0)
        # volume scaling:
        vol = abs(np.linalg.det(newcell))
        newcell *= (self.scaling_volume/vol)**(1./3)

        # Reorder the atoms to follow the atom types in original order
        pos_new = [use_total[n].pop().scaled_position for n in num]

        child = Atoms(numbers=num, scaled_positions=pos_new, 
                     pbc=a1.get_pbc(), cell=newcell)
        return child


    def get_new_individual(self, parents):
        """ The method called by the user that
        returns the paired structure. """
        f, m = parents

        indi = self.cross(f, m)
        desc = 'pairing: {0} {1}'.format(f.info['confid'],
                                         m.info['confid'])
        # It is ok for an operator to return None
        # It means that it could not make a legal offspring
        # within a reasonable amount of time
        if indi is None:
            return indi, desc
        indi = self.initialize_individual(f, indi)
        indi.info['data']['parents'] = [f.info['confid'],
                                        m.info['confid']]
        
        return self.finalize_individual(indi), desc


    def cross(self, a1, a2):
        """Crosses the two atoms objects and returns one"""
        
        if len(a1) != self.n_top: 
            raise ValueError('Wrong size of structure to optimize')
        if len(a1) != len(a2):
            raise ValueError('The two structures do not have the same length')
        
        N = self.n_top

        # Only consider the atoms to optimize
        a1 = a1[len(a1) - N: len(a1)]
        a2 = a2[len(a2) - N: len(a2)]
        
        if not np.array_equal(a1.numbers, a2.numbers):
            err = 'Trying to pair two structures with different stoichiometry'
            raise ValueError(err)

        invalid = True
        counter = 0
        maxcount = 1000
        # Run until a valid pairing is made or 1000 pairings are tested.
        while not_valid and counter < maxcount:

            a1_copy = a1.copy()
            a2_copy = a2.copy()

            # choose one of the 3 lattice vectors:
            direction = randrange(3)

            # shift individuals:
            for a in [a1_copy,a2_copy]:
                cell = a.get_cell()
                for i in range(3):
                    r = random()
                    cond1 = i == direction and r < self.p1
                    cond2 = i != direction and r < self.p2
                    if cond1 or cond2:
                        a.positions += random()*cell[i,:]
                a.wrap()

            # perform the pairing
            fraction = random()
            child = self._get_pairing_(a1_copy, a2_copy, direction=direction,
                                       fraction=fraction)

            # Now checking if the child is a valid candidate

            # Verify whether the atoms are too close or not
            invalid = atoms_too_close(child, self.blmin, 
                                      use_tags=self.use_tags)

            # Verify that the generated structure contains atoms 
            # from both parents
            n1 = -1*np.ones((N,))
            n2 = -1*np.ones((N,))
            scalpos1 = a1_copy.get_scaled_positions()
            scalpos2 = a2_copy.get_scaled_positions()
            scalpos3 = child.get_scaled_positions()
            if self.minfrac is not None:
                nmin1 = int(round(self.minfrac*len(a1_copy)))
                nmin2 = int(round(self.minfrac*len(a2_copy)))
            else:
                nmin1 = 1
                nmin2 = 1
            # Using np.allclose because of some float rounding 
            # when creating the child from the paired cell and positions
            for i in range(N):
                for j in range(N):
                    if np.allclose(scalpos1[j,:],scalpos3[i,:]):
                        n1[i] = j
                        break
                    elif np.allclose(scalpos2[j,:],scalpos3[i,:]):
                        n2[i] = j
                        break
                assert (n1[i] > -1 and n2[i] == -1) or (n1[i] == -1 and
                                                        n2[i] > -1)
            
            if not (len(n1[n1 > -1]) >= nmin1 and len(n2[n2 > -1]) >= nmin2):
                invalid = True

            counter += 1

        if counter == maxcount:
            return None
        
        return child
