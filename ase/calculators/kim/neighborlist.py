from collections import defaultdict

import numpy as np
import kimpy
from kimpy import neighlist
from ase.neighborlist import neighbor_list
from ase import Atom

from .kimpy_wrappers import check_call_wrapper


class NeighborList(object):

    kimpy_arrays = {
        "num_particles": np.intc,
        "coords": np.double,
        "particle_contributing": np.intc,
        "species_code": np.intc,
        "cutoffs": np.double,
        "padding_image_of": np.intc,
        "need_neigh": np.intc,
    }

    def __setattr__(self, name, value):
        if name in self.kimpy_arrays and value is not None:
            value = np.array(value, dtype=self.kimpy_arrays[name])
        self.__dict__[name] = value

    def __init__(
        self,
        neigh_skin_ratio,
        model_influence_dist,
        model_cutoffs,
        padding_not_require_neigh,
        debug,
    ):

        self.skin = neigh_skin_ratio * model_influence_dist
        self.influence_dist = model_influence_dist + self.skin
        self.cutoffs = model_cutoffs + self.skin
        self.padding_need_neigh = not padding_not_require_neigh.all()
        self.debug = debug

        if self.debug:
            print()
            print("Calculator skin: {}".format(self.skin))
            print("Model influence distance:".format(model_influence_dist))
            print(
                "Calculator influence distance (including skin distance): {}"
                "".format(self.influence_dist)
            )
            print("Number of cutoffs: {}".format(model_cutoffs.size))
            print("Model cutoffs: {}".format(model_cutoffs))
            print(
                "Calculator cutoffs (including skin distance): {}"
                "".format(self.cutoffs)
            )
            print(
                "Model needs neighbors of padding atoms: {}"
                "".format(self.padding_need_neigh)
            )
            print()

        # Attributes to be set by subclasses
        self.neigh = None
        self.num_contributing_particles = None
        self.padding_image_of = None
        self.num_particles = None
        self.coords = None
        self.particle_contributing = None
        self.species_code = None
        self.need_neigh = None
        self.last_update_positions = None

    def need_neigh_update(self, atoms, system_changes):
        need_neigh_update = True
        if len(system_changes) == 1 and "positions" in system_changes:
            # only position changes
            if self.last_update_positions is not None:
                a = self.last_update_positions
                b = atoms.positions
                if a.shape == b.shape:
                    delta = np.linalg.norm(a - b, axis=1)
                    # indices of the two largest elements
                    ind = np.argpartition(delta, -2)[-2:]
                    if sum(delta[ind]) <= self.skin:
                        need_neigh_update = False

        return need_neigh_update

    def clean(self):
        pass


class ASENeighborList(NeighborList):
    def __init__(
        self,
        compute_args,
        neigh_skin_ratio,
        model_influence_dist,
        model_cutoffs,
        padding_not_require_neigh,
        debug,
    ):
        super().__init__(
            neigh_skin_ratio,
            model_influence_dist,
            model_cutoffs,
            padding_not_require_neigh,
            debug,
        )

        self.neigh = {}
        compute_args.set_callback(
            kimpy.compute_callback_name.GetNeighborList, self.get_neigh, self.neigh
        )

    @staticmethod
    def get_neigh(data, cutoffs, neighbor_list_index, particle_number):
        """Retrieves the neighbors of each atom using ASE's native neighbor
        list library
        """
        # We can only return neighbors of particles that were stored
        number_of_particles = data["num_particles"]
        if particle_number >= number_of_particles or particle_number < 0:
            return (np.array([]), 1)

        neighbors = data["neighbors"][neighbor_list_index][particle_number]
        return (neighbors, 0)

    def build(self, atoms):
        """Build the ASE neighbor list and return an Atoms object with
        all of the neighbors added.  First a neighbor list is created
        from ase.neighbor_list, having only information about the
        neighbors of the original atoms.  If neighbors of padding atoms
        are required, they are calculated using information from the
        first neighbor list.
        """
        syms = atoms.get_chemical_symbols()
        num_atoms = len(atoms)
        i, j, D, S, dists = neighbor_list("ijDSd", atoms, self.influence_dist)

        # Get coordinates for all neighbors (this has overlapping positions)
        A = atoms.get_positions()[i] + D

        # Make the neighbor list ready for KIM
        ac = atoms.copy()
        used = dict()

        # Variables below only include information for the neighbors (padding)
        padding_image_of = []
        neighbor_shifts = []

        neigh_dict = defaultdict(list)
        neigh_dists = defaultdict(list)

        # Loop over all neighbor pairs
        for k in range(len(i)):
            shift_tuple = tuple(S[k])
            t = (j[k],) + shift_tuple
            if shift_tuple == (0, 0, 0):
                # In unit cell
                neigh_dict[i[k]].append(j[k])
                neigh_dists[i[k]].append(dists[k])
                if t not in used:
                    used[t] = j[k]
            else:
                # Not in unit cell
                if t not in used:
                    # Add the neighbor as a padding atom
                    used[t] = len(ac)
                    ac.append(Atom(syms[j[k]], position=A[k]))
                    padding_image_of.append(j[k])
                    neighbor_shifts.append(S[k])
                neigh_dict[i[k]].append(used[t])
                neigh_dists[i[k]].append(dists[k])
        neighbor_list_size = num_atoms

        # Add neighbors of padding atoms if the potential requires them
        if self.padding_need_neigh:
            neighbor_list_size = len(ac)
            inv_used = dict((v, k) for k, v in used.items())
            # Loop over all the neighbors (k)
            # and the image of that neighbor in the cell (neigh)
            for k, neigh in enumerate(padding_image_of):
                # Shift from original atom in cell to neighbor
                shift = neighbor_shifts[k]
                for org_neigh, org_dist in zip(neigh_dict[neigh], neigh_dists[neigh]):
                    # Get the shift of the neighbor of the original atom
                    org_shift = inv_used[org_neigh][1:]

                    # Apply sum of original shift and current shift
                    # to neighbors of original atom
                    tot_shift = org_shift + shift

                    # Get the image in the cell of the original neighbor
                    if org_neigh <= num_atoms - 1:
                        org_neigh_image = org_neigh
                    else:
                        org_neigh_image = padding_image_of[org_neigh - num_atoms]

                    # If the original image with the total shift has been
                    # used before then it is also a neighbor of this atom
                    tt = (org_neigh_image,) + tuple(tot_shift)
                    if tt in used:
                        neigh_dict[k + num_atoms].append(used[tt])
                        neigh_dists[k + num_atoms].append(org_dist)

        neigh_lists = []
        for cut in self.cutoffs:
            neigh_list = [
                np.array(neigh_dict[k], dtype=np.intc)[neigh_dists[k] <= cut]
                for k in range(neighbor_list_size)
            ]
            neigh_lists.append(neigh_list)

        self.padding_image_of = padding_image_of

        self.neigh["neighbors"] = neigh_lists
        self.neigh["num_particles"] = neighbor_list_size

        return ac

    def update(self, atoms, species_map):
        """Create the neighbor list along with the other required
        parameters (which are stored as instance attributes). The
        required parameters are:

            - num_particles
            - coords
            - particle_contributing
            - species_code

        Note that the KIM API requires a neighbor list that has indices
        corresponding to each atom.
        """

        # Information of original atoms
        self.num_contributing_particles = len(atoms)

        ac = self.build(atoms)

        # Save the number of atoms and all their neighbors and positions
        N = len(ac)
        num_padding = N - self.num_contributing_particles
        self.num_particles = [N]
        self.coords = ac.get_positions()

        # Save which coordinates are from original atoms and which are from
        # neighbors using a mask
        indices_mask = [1] * self.num_contributing_particles + [0] * num_padding
        self.particle_contributing = indices_mask

        # species support and code
        try:
            self.species_code = [species_map[s] for s in ac.get_chemical_symbols()]
        except KeyError as e:
            raise RuntimeError("Species not supported by KIM model; {}".format(str(e)))

        self.last_update_positions = atoms.get_positions()

        if self.debug:
            print("Debug: called update_ase_neigh")
            print()


class KimpyNeighborList(NeighborList):
    def __init__(
        self,
        compute_args,
        neigh_skin_ratio,
        model_influence_dist,
        model_cutoffs,
        padding_not_require_neigh,
        debug,
    ):
        super().__init__(
            neigh_skin_ratio,
            model_influence_dist,
            model_cutoffs,
            padding_not_require_neigh,
            debug,
        )

        self.neigh = neighlist.initialize()
        compute_args.set_callback_pointer(
            kimpy.compute_callback_name.GetNeighborList,
            neighlist.get_neigh_kim(),
            self.neigh,
        )

    @check_call_wrapper
    def build(self):
        return neighlist.build(
            self.neigh, self.coords, self.influence_dist, self.cutoffs, self.need_neigh
        )

    @check_call_wrapper
    def create_paddings(
        self, cell, pbc, contributing_coords, contributing_species_code
    ):
        # Cast things passed through kimpy to numpy arrays
        cell = np.asarray(cell, dtype=np.double)
        pbc = np.asarray(pbc, dtype=np.intc)
        contributing_coords = np.asarray(contributing_coords, dtype=np.double)

        return neighlist.create_paddings(
            self.influence_dist,
            cell,
            pbc,
            contributing_coords,
            contributing_species_code,
        )

    def update(self, atoms, species_map):
        """Create the neighbor list along with the other required
        parameters (which are stored as instance attributes). The
        required parameters are:

            - num_particles
            - coords
            - particle_contributing
            - species_code

        Note that the KIM API requires a neighbor list that has indices
        corresponding to each atom.
        """

        # get info from Atoms object
        cell = np.asarray(atoms.get_cell(), dtype=np.double)
        pbc = np.asarray(atoms.get_pbc(), dtype=np.intc)
        contributing_coords = np.asarray(atoms.get_positions(), dtype=np.double)
        self.num_contributing_particles = atoms.get_global_number_of_atoms()
        num_contributing = self.num_contributing_particles

        # species support and code
        try:
            contributing_species_code = np.array(
                [species_map[s] for s in atoms.get_chemical_symbols()], dtype=np.intc
            )
        except KeyError as e:
            raise RuntimeError("Species not supported by KIM model; {}".format(str(e)))

        if pbc.any():  # need padding atoms
            # create padding atoms

            padding_coords, padding_species_code, self.padding_image_of = self.create_paddings(
                cell, pbc, contributing_coords, contributing_species_code
            )
            num_padding = padding_species_code.size

            self.num_particles = [num_contributing + num_padding]
            self.coords = np.concatenate((contributing_coords, padding_coords))
            self.species_code = np.concatenate(
                (contributing_species_code, padding_species_code)
            )
            self.particle_contributing = [1] * num_contributing + [0] * num_padding
            self.need_neigh = [1] * self.num_particles[0]
            if not self.padding_need_neigh:
                self.need_neigh[num_contributing:] = 0

        else:  # do not need padding atoms
            self.padding_image_of = []
            self.num_particles = [num_contributing]
            self.coords = contributing_coords
            self.species_code = contributing_species_code
            self.particle_contributing = [1] * num_contributing
            self.need_neigh = self.particle_contributing

        # create neighborlist
        self.build()

        self.last_update_positions = atoms.get_positions()

        if self.debug:
            print("Debug: called update_kimpy_neigh")
            print()

    def clean(self):
        neighlist.clean(self.neigh)
