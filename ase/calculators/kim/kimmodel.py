"""
ASE Calculator for interatomic models compatible with the Knowledgebase
of Interatomic Models (KIM) application programming interface (API).
Written by:

Mingjian Wen
Daniel S. Karls
University of Minnesota
"""
import numpy as np

from ase.calculators.calculator import Calculator
from ase.calculators.calculator import equal

from . import kimpy_wrappers
from . import neighborlist


class KIMModelData(object):
    """Initializes and subsequently stores the KIM API Model object, KIM
    API ComputeArguments object, and the neighbor list object used by
    instances of KIMModelCalculator.  Also stores the arrays which are
    registered in the KIM API and which are used to communicate with the
    model.
    """

    def __init__(self, model_name, ase_neigh, neigh_skin_ratio, debug=False):
        self.model_name = model_name
        self.ase_neigh = ase_neigh
        self.debug = debug

        self.kim_initialized = False
        self.neigh_initialized = False

        # model object, compute arguments object, neighbor list object
        self.kim_model = None
        self.compute_args = None
        self.neigh = None

        # initialize KIM API Portable Model object and ComputeArguments object
        self.init_kim()

        # set cutoff
        model_influence_dist = self.kim_model.get_influence_distance()
        model_cutoffs, padding_not_require_neigh = (
            self.kim_model.get_neighbor_list_cutoffs_and_hints()
        )

        self.species_map = self.create_species_map()

        # initialize neighbor list object
        self.init_neigh(
            neigh_skin_ratio,
            model_influence_dist,
            model_cutoffs,
            padding_not_require_neigh,
        )

    def init_kim(self):
        """Create the KIM API Model object and KIM API ComputeArguments
        object
        """

        if self.kim_initialized:
            return

        self.kim_model = kimpy_wrappers.PortableModel(self.model_name, self.debug)

        # KIM API model object is what actually creates/destroys the ComputeArguments
        # object, so we must pass it as a parameter
        self.compute_args = self.kim_model.compute_arguments_create()

        self.kim_initialized = True

    def init_neigh(
        self,
        neigh_skin_ratio,
        model_influence_dist,
        model_cutoffs,
        padding_not_require_neigh,
    ):

        """Initialize neighbor list, either an ASE-native neighborlist
        or one created using the neighlist module in kimpy
        """
        if self.ase_neigh:
            self.neigh = neighborlist.ASENeighborList(
                self.compute_args,
                neigh_skin_ratio,
                model_influence_dist,
                model_cutoffs,
                padding_not_require_neigh,
                self.debug,
            )
        else:
            self.neigh = neighborlist.KimpyNeighborList(
                self.compute_args,
                neigh_skin_ratio,
                model_influence_dist,
                model_cutoffs,
                padding_not_require_neigh,
                self.debug,
            )

        self.neigh_initialized = True

    def update_kim_coords(self, atoms):
        """Update atomic positions in self.coords, which where the KIM
        API will look to find them in order to pass them to the model.
        """
        if self.padding_image_of.size != 0:
            disp_contrib = atoms.positions - self.coords[: len(atoms)]
            disp_pad = disp_contrib[self.padding_image_of]
            self.coords += np.concatenate((disp_contrib, disp_pad))
        else:
            np.copyto(self.coords, atoms.positions)

        if self.debug:
            print("Debug: called update_kim_coords")
            print()

    def update_compute_args_pointers(self, energy, forces):
        self.compute_args.update(
            self.num_particles,
            self.species_code,
            self.particle_contributing,
            self.coords,
            energy,
            forces,
        )

    def create_species_map(self):
        """Get all the supported species of the KIM model and the
        corresponding integer codes used by the model

        Returns
        -------
        species_map : dict
            key : str
                chemical symbols (e.g. "Ar")
            value : int
                species integer code (e.g. 1)
        """
        supported_species, codes = self.get_model_supported_species_and_codes()
        species_map = dict()
        for i, s in enumerate(supported_species):
            species_map[s] = codes[i]
            if self.debug:
                print("Species {} is supported and its code is: {}".format(s, codes[i]))

        return species_map

    def clean_neigh(self):
        """If the neighbor list method being used is the one in the
        kimpy neighlist module, deallocate its memory
        """
        if self.neigh_initialized:
            self.neigh.clean()
            self.neigh_initialized = False

    def clean_kim(self):
        """Deallocate the memory allocated to the KIM API Model object
        and KIM API ComputeArguments object
        """
        if self.kim_initialized:
            self.kim_model.compute_arguments_destroy(self.compute_args)
            self.kim_model.destroy()
            self.kim_initialized = False

    def clean(self):
        """Deallocate the KIM API Model object, KIM API ComputeArguments
        object, and, if applicable, the neighbor list object
        """
        self.clean_neigh()
        self.clean_kim()

    def __del__(self):
        self.clean()

    @property
    def padding_image_of(self):
        return self.neigh.padding_image_of

    @property
    def num_particles(self):
        return self.neigh.num_particles

    @property
    def coords(self):
        return self.neigh.coords

    @property
    def particle_contributing(self):
        return self.neigh.particle_contributing

    @property
    def species_code(self):
        return self.neigh.species_code

    @property
    def get_model_supported_species_and_codes(self):
        return self.kim_model.get_model_supported_species_and_codes


class KIMModelCalculator(Calculator):
    """Calculator that works with KIM Portable Models (PMs).

    Calculator that carries out direct communication between ASE and a
    KIM Portable Model (PM) through the kimpy library (which provides a
    set of python bindings to the KIM API).

    Parameters
    ----------
    model_name : str
      The unique identifier assigned to the interatomic model (for
      details, see https://openkim.org/doc/schema/kim-ids)

    ase_neigh : bool, optional
      False (default): Use kimpy's neighbor list library

      True: Use ASE's internal neighbor list mechanism (usually slower
      than the kimpy neighlist library)

    neigh_skin_ratio : float, optional
      Used to determine the neighbor list cutoff distance, r_neigh,
      through the relation r_neigh = (1 + neigh_skin_ratio) * rcut,
      where rcut is the model's influence distance. (Default: 0.2)

    release_GIL : bool, optional
      Whether to release python GIL.  Releasing the GIL allows a KIM
      model to run with multiple concurrent threads. (Default: False)

    debug : bool, optional
      If True, detailed information is printed to stdout. (Default:
      False)
    """

    implemented_properties = ["energy", "forces", "stress"]

    def __init__(
        self,
        model_name,
        ase_neigh=False,
        neigh_skin_ratio=0.2,
        release_GIL=False,
        debug=False,
        *args,
        **kwargs
    ):
        super().__init__(*args, **kwargs)

        self.model_name = model_name
        self.release_GIL = release_GIL
        self.debug = debug

        # neigh attributes
        if neigh_skin_ratio < 0:
            raise ValueError('Argument "neigh_skin_ratio" must be non-negative')

        # model output
        self.energy = None
        self.forces = None

        # create KIMModelData object. This will take care of creating and storing the
        # KIM API Model object, KIM API ComputeArguments object, and the neighbor list
        # object that our calculator needs
        self.kimmodeldata = KIMModelData(
            self.model_name, ase_neigh, neigh_skin_ratio, self.debug
        )

    def __enter__(self):
        return self

    def __exit__(self, exc_type, value, traceback):
        if exc_type is None:
            if self.kimmodeldata is not None:
                # Explicitly deallocate all three objects held by the KIMModelData
                # instance referenced by our calculator
                self.kimmodeldata.clean()
        else:
            return False  # reraise exception

    def __repr__(self):
        return "KIMModelCalculator(model_name={})".format(self.model_name)

    def calculate(
        self,
        atoms=None,
        properties=["energy", "forces", "stress"],
        system_changes=["positions", "numbers", "cell", "pbc"],
    ):
        """
        Inherited method from the ase Calculator class that is called by
        get_property()

        Parameters
        ----------
        atoms : Atoms
            Atoms object whose properties are desired

        properties : list of str
            List of what needs to be calculated.  Can be any combination
            of 'energy', 'forces' and 'stress'.

        system_changes : list of str
            List of what has changed since last calculation.  Can be any
            combination of these six: 'positions', 'numbers', 'cell',
            and 'pbc'.
        """

        Calculator.calculate(self, atoms, properties, system_changes)

        # update KIM API input data and neighbor list if necessary
        if system_changes:
            if self.need_neigh_update(atoms, system_changes):
                self.update_neigh(atoms, self.species_map)
                self.energy = np.array([0.0], dtype=np.double)
                self.forces = np.zeros([self.num_particles[0], 3], dtype=np.double)
                self.update_compute_args_pointers(self.energy, self.forces)
            else:
                self.update_kim_coords(atoms)

            self.kim_model.compute(self.compute_args, self.release_GIL)

        energy = self.energy[0]
        forces = self.assemble_padding_forces()

        try:
            volume = atoms.get_volume()
            stress = self.compute_virial_stress(self.forces, self.coords, volume)
        except ValueError:  # volume cannot be computed
            stress = None

        # return values
        self.results["energy"] = energy
        self.results["forces"] = forces
        self.results["stress"] = stress

    def check_state(self, atoms, tol=1e-15):
        return self.compare_atoms(self.atoms, atoms)

    @staticmethod
    def compare_atoms(atoms1, atoms2, tol=1e-15):
        """Check for system changes since last calculation. Note that
        this is an override of Calculator.compare_atoms and differs in
        that the magnetic moments and charges are not checked because
        the KIM API does not (currently) support these.
        """
        if atoms1 is None:
            return ["positions", "numbers", "cell", "pbc"]
        else:
            system_changes = []
            if not equal(atoms1.positions, atoms2.positions, tol):
                system_changes.append("positions")
            if not equal(atoms1.numbers, atoms2.numbers):
                system_changes.append("numbers")
            if not equal(atoms1.cell, atoms2.cell, tol):
                system_changes.append("cell")
            if not equal(atoms1.pbc, atoms2.pbc):
                system_changes.append("pbc")

        return system_changes

    def assemble_padding_forces(self):
        """
        Assemble forces on padding atoms back to contributing atoms.

        Parameters
        ----------
        forces : 2D array of doubles
            Forces on both contributing and padding atoms

        num_contrib:  int
            Number of contributing atoms

        padding_image_of : 1D array of int
            Atom number, of which the padding atom is an image


        Returns
        -------
            Total forces on contributing atoms.
        """

        total_forces = np.array(self.forces[: self.num_contributing_particles])

        if self.padding_image_of.size != 0:
            pad_forces = self.forces[self.num_contributing_particles :]
            for f, org_index in zip(pad_forces, self.padding_image_of):
                total_forces[org_index] += f

        return total_forces

    @staticmethod
    def compute_virial_stress(forces, coords, volume):
        """Compute the virial stress in Voigt notation.

        Parameters
        ----------
        forces : 2D array
            Partial forces on all atoms (padding included)

        coords : 2D array
            Coordinates of all atoms (padding included)

        volume : float
            Volume of cell

        Returns
        -------
        stress : 1D array
            stress in Voigt order (xx, yy, zz, yz, xz, xy)
        """
        stress = np.zeros(6)
        stress[0] = -np.dot(forces[:, 0], coords[:, 0]) / volume
        stress[1] = -np.dot(forces[:, 1], coords[:, 1]) / volume
        stress[2] = -np.dot(forces[:, 2], coords[:, 2]) / volume
        stress[3] = -np.dot(forces[:, 1], coords[:, 2]) / volume
        stress[4] = -np.dot(forces[:, 0], coords[:, 2]) / volume
        stress[5] = -np.dot(forces[:, 0], coords[:, 1]) / volume

        return stress

    @property
    def update_compute_args_pointers(self):
        return self.kimmodeldata.update_compute_args_pointers

    @property
    def get_model_supported_species_and_codes(self):
        return self.kimmodeldata.get_model_supported_species_and_codes

    @property
    def kim_model(self):
        return self.kimmodeldata.kim_model

    @property
    def compute_args(self):
        return self.kimmodeldata.compute_args

    @property
    def num_particles(self):
        return self.kimmodeldata.num_particles

    @property
    def coords(self):
        return self.kimmodeldata.coords

    @property
    def padding_image_of(self):
        return self.kimmodeldata.padding_image_of

    @property
    def update_kim_coords(self):
        return self.kimmodeldata.update_kim_coords

    @property
    def species_map(self):
        return self.kimmodeldata.species_map

    @property
    def neigh(self):
        return self.kimmodeldata.neigh

    @property
    def num_contributing_particles(self):
        return self.neigh.num_contributing_particles

    @property
    def need_neigh_update(self):
        return self.neigh.need_neigh_update

    @property
    def update_neigh(self):
        return self.neigh.update
