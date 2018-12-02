""" A collection of mutations that can be used. """

import numpy as np
from random import random, randrange
from math import ceil, cos, sin, pi
from ase.ga.utilities import atoms_too_close
from ase.ga.utilities import atoms_too_close_two_sets
from ase.ga.offspring_creator import OffspringCreator
from ase import Atoms


class RattleMutation(OffspringCreator):
    """ An implementation of the rattle mutation as described in
        R.L. Johnston Dalton Transactions, Vol. 22,
        No. 22. (2003), pp. 4193-4207

        Parameters:

        blmin: Dictionary defining the minimum distance between atoms
        after the rattle.
        n_top: Number of atoms optimized by the GA.
        rattle_strength: Strength with which the atoms are moved.
        rattle_prop: The probability with which each atom is rattled.
        test_dist_to_slab: whether to also make sure that the distances
                  between the atoms and the slab satisfy the blmin.
        use_tags: if True, the atomic tags will be used to preserve
                  molecular identity. Same-tag atoms will then be
                  displaced collectively, so that the internal
                  geometry is preserved.
    """
    def __init__(self, blmin, n_top, rattle_strength=0.8,
                 rattle_prop=0.4, test_dist_to_slab=True, use_tags=False,
                 verbose=False):
        OffspringCreator.__init__(self, verbose)
        self.blmin = blmin
        self.n_top = n_top
        self.rattle_strength = rattle_strength
        self.rattle_prop = rattle_prop
        self.test_dist_to_slab = test_dist_to_slab
        self.use_tags = use_tags
        self.descriptor = 'RattleMutation'
        self.min_inputs = 1

    def get_new_individual(self, parents):
        f = parents[0]

        indi = self.mutate(f)
        if indi is None:
            return indi, 'mutation: rattle'
            
        indi = self.initialize_individual(f, indi)
        indi.info['data']['parents'] = [f.info['confid']]

        return self.finalize_individual(indi), 'mutation: rattle'

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
            for tag in np.unique(tags):
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

        
class PermutationMutation(OffspringCreator):
    """Mutation that permutes a percentage of the atom types in the cluster.

       Parameters:

       n_top: Number of atoms optimized by the GA.
       probability: The probability with which an atom is permuted.
    """

    def __init__(self, n_top, probability=0.33, verbose=False):
        OffspringCreator.__init__(self, verbose)
        self.n_top = n_top
        self.probability = probability
        self.descriptor = 'PermutationMutation'
        self.min_inputs = 1

    def get_new_individual(self, parents):
        f = parents[0]

        indi = self.mutate(f)
        if indi is None:
            return indi, 'mutation: permutation'
            
        indi = self.initialize_individual(f, indi)
        indi.info['data']['parents'] = [f.info['confid']]

        return self.finalize_individual(indi), 'mutation: permutation'

    def mutate(self, atoms):
        """ Does the actual mutation. """
        a = atoms.copy()
        s = self.n_top
        p = a.get_positions()[-s:]
        n = a.numbers[-s:]
        n_un = list(set(n))
        assert len(n_un) > 1, 'Permutations with one atomic type is not valid'
        m = int(ceil(float(s) * self.probability / 2.))
        for _ in range(m):
            i = j = 0
            while n[i] == n[j]:
                i = randrange(0, s)
                j = randrange(0, s)
            t = p[i].copy()
            p[i] = p[j].copy()
            p[j] = t
        p_tot = a.get_positions()
        p_tot[-s:] = p
        a.set_positions(p_tot)
        return a


class MirrorMutation(OffspringCreator):
    """ A mirror mutation, as described in
        TO BE PUBLISHED.
        This mutation mirrors half of the cluster in a
        randomly oriented cutting plane discarding the other half.

        Parameters:
        blmin: Dictionary defining the minimum allowed
        distance between atoms.
        n_top: Number of atoms the GA optimizes.
        reflect: Defines if the mirrored half is also reflected
        perpendicular to the mirroring plane.

    """
    def __init__(self, blmin, n_top, reflect=False, verbose=False):
        OffspringCreator.__init__(self, verbose)
        self.blmin = blmin
        self.n_top = n_top
        self.reflect = reflect
        self.descriptor = 'MirrorMutation'
        self.min_inputs = 1

    def get_new_individual(self, parents):
        f = parents[0]

        indi = self.mutate(f)
        if indi is None:
            return indi, 'mutation: mirror'
            
        indi = self.initialize_individual(f, indi)
        indi.info['data']['parents'] = [f.info['confid']]

        return self.finalize_individual(indi), 'mutation: mirror'

    def mutate(self, atoms):
        """ Do the mutation of the atoms input. """

        reflect = self.reflect
        tc = True
        slab = atoms[0:len(atoms) - self.n_top]
        top = atoms[len(atoms) - self.n_top: len(atoms)]
        num = top.numbers
        unique_types = list(set(num))
        nu = dict()
        for u in unique_types:
            nu[u] = sum(num == u)
            
        n_tries = 1000
        counter = 0
        changed = False

        while tc and counter < n_tries:
            counter += 1
            cand = top.copy()
            pos = cand.get_positions()

            cm = np.average(top.get_positions(), axis=0)

            # first select a randomly oriented cutting plane
            theta = pi * random()
            phi = 2. * pi * random()
            n = (cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta))
            n = np.array(n)

            # Calculate all atoms signed distance to the cutting plane
            D = []
            for (i, p) in enumerate(pos):
                d = np.dot(p - cm, n)
                D.append((i, d))

            # Sort the atoms by their signed distance
            D.sort(key=lambda x: x[1])
            nu_taken = dict()

            # Select half of the atoms needed for a full cluster
            p_use = []
            n_use = []
            for (i, d) in D:
                if num[i] not in nu_taken.keys():
                    nu_taken[num[i]] = 0
                if nu_taken[num[i]] < nu[num[i]] / 2.:
                    p_use.append(pos[i])
                    n_use.append(num[i])
                    nu_taken[num[i]] += 1

            # calculate the mirrored position and add these.
            pn = []
            for p in p_use:
                pt = p - 2. * np.dot(p - cm, n) * n
                if reflect:
                    pt = -pt + 2 * cm + 2 * n * np.dot(pt - cm, n)
                pn.append(pt)

            n_use.extend(n_use)
            p_use.extend(pn)

            # In the case of an uneven number of
            # atoms we need to add one extra
            for n in nu.keys():
                if nu[n] % 2 == 0:
                    continue
                while sum(n_use == n) > nu[n]:
                    for i in range(int(len(n_use) / 2), len(n_use)):
                        if n_use[i] == n:
                            del p_use[i]
                            del n_use[i]
                            break
                assert sum(n_use == n) == nu[n]

            # Make sure we have the correct number of atoms
            # and rearrange the atoms so they are in the right order
            for i in range(len(n_use)):
                if num[i] == n_use[i]:
                    continue
                for j in range(i + 1, len(n_use)):
                    if n_use[j] == num[i]:
                        tn = n_use[i]
                        tp = p_use[i]
                        n_use[i] = n_use[j]
                        p_use[i] = p_use[j]
                        p_use[j] = tp
                        n_use[j] = tn

            # Finally we check that nothing is too close in the end product.
            cand = Atoms(num, p_use, cell=slab.get_cell(), pbc=slab.get_pbc())
            tc = atoms_too_close(cand, self.blmin)
            if tc:
                continue
            tc = atoms_too_close_two_sets(slab, cand, self.blmin)
            if not changed and counter > n_tries // 2:
                reflect = not reflect
                changed = True
            tot = slab + cand
        if counter == n_tries:
            return None
        return tot
