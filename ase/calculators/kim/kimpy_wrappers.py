"""
Wrappers that provide a minimal interface to kimpy methods and objects

Daniel S. Karls
University of Minnesota
"""

import functools

import kimpy
from .kim import KIMCalculatorError

def check_call(f, *args):
    """
    Wrapper for functions that checks error codes, since many of the
    functions calls are actually python bindings to functions written in
    other languages.  Functions are assumed to return either an integer
    error code or a tuple whose last element is an integer error code.
    """

    def _check_error(error, msg):
        if error != 0 and error is not None:
            raise KIMCalculatorError('Calling "{}" failed.'.format(msg))

    ret = f(*args)

    if isinstance(ret, int):
        # Only an error code was returned
        _check_error(ret, f.__name__)
    else:
        # An error code plus other variables were returned
        error = ret[-1]
        _check_error(error, f.__name__)

        if len(ret[:-1]) == 1:
            # Pick the single remaining element out of the tuple
            return ret[0]
        else:
            # Return the tuple containing the rest of the elements
            return ret[:-1]


def check_call_wrapper(func):
    @functools.wraps(func)
    def myfunc(*args, **kwargs):
        return check_call(func, *args)

    return myfunc


# kimpy methods wrapped in ``check_error``
collections_create = functools.partial(check_call, kimpy.collections.create)
model_create = functools.partial(check_call, kimpy.model.create)
simulator_model_create = functools.partial(check_call, kimpy.simulator_model.create)
get_species_name = functools.partial(check_call, kimpy.species_name.get_species_name)

# kimpy attributes (here to avoid importing kimpy in higher-level modules)
collection_item_type_portableModel = kimpy.collection_item_type.portableModel


class ModelCollections(object):
    def __init__(self):
        self.initialized = False
        self.collection = collections_create()
        self.initialized = True

    def __del__(self):
        self.destroy()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, value, traceback):
        if exc_type is None:
            self.destroy()
        else:
            return False  # reraise exception

    def get_item_type(self, model_name):
        try:
            model_type = check_call(self.collection.get_item_type, model_name)
        except KIMCalculatorError as e:
            msg = (
                "ERROR: Could not find model {} installed in any of the KIM API model "
                "collections on this system.  See "
                "https://openkim.org/doc/usage/obtaining-models/ for instructions on "
                "installing models.".format(model_name)
            )
            raise e(msg)

        return model_type

    def destroy(self):
        if self.initialized:
            kimpy.collections.destroy(self.collection)
            self.initialized = False


class PortableModel(object):
    """ Provides a minimal interface (only what we require here) to the
    KIM API Model object that we instantiate
    """

    def __init__(self, model_name, debug):
        self.model_name = model_name
        self.debug = debug

        # create KIM API Model object
        units_accepted, self.kim_model = model_create(
            kimpy.numbering.zeroBased,
            kimpy.length_unit.A,
            kimpy.energy_unit.eV,
            kimpy.charge_unit.e,
            kimpy.temperature_unit.K,
            kimpy.time_unit.ps,
            self.model_name,
        )

        if not units_accepted:
            raise KIMCalculatorError("Requested units not accepted in kimpy.model.create")

        if self.debug:
            l_unit, e_unit, c_unit, te_unit, ti_unit = self.kim_model.get_units()
            print("Length unit is: {}".format(l_unit))
            print("Energy unit is: {}".format(e_unit))
            print("Charge unit is: {}".format(c_unit))
            print("Temperature unit is: {}".format(te_unit))
            print("Time unit is: {}".format(ti_unit))
            print()

    def get_model_supported_species_and_codes(self):
        """Get all the supported species and corresponding integer codes
        for the KIM Portable Model.

        Returns
        -------
        species : list of str
            Abbreviated chemical symbols of all species the mmodel
            supports (e.g. ["Mo", "S"])

        codes : list of int
            Integer codes used by the model for each species (order
            corresponds to the order of ``species``)
        """
        species = []
        codes = []
        num_kim_species = kimpy.species_name.get_number_of_species_names()

        for i in range(num_kim_species):
            species_name = get_species_name(i)
            species_support, code = self.get_species_support_and_code(species_name)

            if species_support:
                species.append(str(species_name))
                codes.append(code)

        return species, codes

    @check_call_wrapper
    def compute(self, compute_args_wrapped, release_GIL):
        return self.kim_model.compute(compute_args_wrapped.compute_args, release_GIL)

    @check_call_wrapper
    def get_species_support_and_code(self, species_name):
        return self.kim_model.get_species_support_and_code(species_name)

    def get_influence_distance(self):
        return self.kim_model.get_influence_distance()

    def get_neighbor_list_cutoffs_and_hints(self):
        return self.kim_model.get_neighbor_list_cutoffs_and_hints()

    def compute_arguments_create(self):
        return ComputeArguments(self, self.debug)

    def compute_arguments_destroy(self, compute_args_wrapped):
        compute_args_wrapped.destroy()

    def destroy(self):
        kimpy.model.destroy(self.kim_model)


class ComputeArguments(object):
    """
    KIM API ComputeArguments object configured for ASE
    """

    def __init__(self, kim_model_wrapped, debug):
        self.kim_model_wrapped = kim_model_wrapped
        self.debug = debug

        # create KIM API ComputeArguments object
        self.compute_args = check_call(
            self.kim_model_wrapped.kim_model.compute_arguments_create
        )

        # check compute arguments
        kimpy_arg_name = kimpy.compute_argument_name
        num_arguments = kimpy_arg_name.get_number_of_compute_argument_names()
        if self.debug:
            print("Number of compute_args: {}".format(num_arguments))

        for i in range(num_arguments):
            name = check_call(kimpy_arg_name.get_compute_argument_name, i)
            dtype = check_call(kimpy_arg_name.get_compute_argument_data_type, name)

            arg_support = self.get_argument_support_status(name)

            if self.debug:
                print(
                    "Compute Argument name {:21} is of type {:7} and has support "
                    "status {}".format(*[str(x) for x in [name, dtype, arg_support]])
                )

            # the simulator can handle energy and force from a kim model
            # virial is computed within the calculator
            if arg_support == kimpy.support_status.required:
                if (
                    name != kimpy.compute_argument_name.partialEnergy
                    and name != kimpy.compute_argument_name.partialForces
                ):
                    raise KIMCalculatorError("Unsupported required ComputeArgument {}".format(name))

        # check compute callbacks
        callback_name = kimpy.compute_callback_name
        num_callbacks = callback_name.get_number_of_compute_callback_names()
        if self.debug:
            print()
            print("Number of callbacks: {}".format(num_callbacks))

        for i in range(num_callbacks):
            name = check_call(callback_name.get_compute_callback_name, i)

            support_status = self.get_callback_support_status(name)

            if self.debug:
                print(
                    "Compute callback {:17} has support status {}".format(
                        str(name), support_status
                    )
                )

            # cannot handle any "required" callbacks
            if support_status == kimpy.support_status.required:
                raise KIMCalculatorError("Unsupported required ComputeCallback: {}".format(name))

    @check_call_wrapper
    def set_argument_pointer(self, compute_arg_name, data_object):
        return self.compute_args.set_argument_pointer(compute_arg_name, data_object)

    @check_call_wrapper
    def get_argument_support_status(self, name):
        return self.compute_args.get_argument_support_status(name)

    @check_call_wrapper
    def get_callback_support_status(self, name):
        return self.compute_args.get_callback_support_status(name)

    @check_call_wrapper
    def set_callback(self, compute_callback_name, callback_function, data_object):
        return self.compute_args.set_callback(
            compute_callback_name, callback_function, data_object
        )

    @check_call_wrapper
    def set_callback_pointer(self, compute_callback_name, callback, data_object):
        return self.compute_args.set_callback_pointer(
            compute_callback_name, callback, data_object
        )

    def update(
        self, num_particles, species_code, particle_contributing, coords, energy, forces
    ):
        """ Register model input and output in the kim_model object."""
        compute_arg_name = kimpy.compute_argument_name
        set_argument_pointer = self.set_argument_pointer

        set_argument_pointer(compute_arg_name.numberOfParticles, num_particles)
        set_argument_pointer(compute_arg_name.particleSpeciesCodes, species_code)
        set_argument_pointer(
            compute_arg_name.particleContributing, particle_contributing
        )
        set_argument_pointer(compute_arg_name.coordinates, coords)
        set_argument_pointer(compute_arg_name.partialEnergy, energy)
        set_argument_pointer(compute_arg_name.partialForces, forces)

        if self.debug:
            print("Debug: called update_kim")
            print()

    @check_call_wrapper
    def destroy(self):
        return self.kim_model_wrapped.kim_model.compute_arguments_destroy(
            self.compute_args
        )


class SimulatorModel(object):
    def __init__(self, model_name):
        self.initialized = False

        # Create a KIM API simulator Model object for this model
        self.model_name = model_name
        self.simulator_model = simulator_model_create(self.model_name)

        # Need to close template map to access simulator model metadata
        self.simulator_model.close_template_map()

        self.initialized = True

    def __del__(self):
        self.destroy()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, value, traceback):
        if exc_type is None:
            self.destroy()
        else:
            return False  # reraise exception

    @property
    def simulator_name(self):
        simulator_name, _ = self.simulator_model.get_simulator_name_and_version()
        return simulator_name

    @property
    def num_supported_species(self):
        num_supported_species = self.simulator_model.get_number_of_supported_species()
        if num_supported_species == 0:
            raise KIMCalculatorError(
                "ERROR: Unable to determine supported species of "
                "simulator model {}.".format(self.model_name)
            )
        else:
            return num_supported_species

    @property
    def supported_species(self):
        supported_species = []
        for spec_code in range(self.num_supported_species):
            species = check_call(self.simulator_model.get_supported_species, spec_code)
            supported_species.append(species)

        return tuple(supported_species)

    @property
    def num_metadata_fields(self):
        return self.simulator_model.get_number_of_simulator_fields()

    @property
    def metadata(self):
        sm_metadata_fields = {}
        for field in range(self.num_metadata_fields):
            extent, field_name = check_call(
                self.simulator_model.get_simulator_field_metadata, field
            )
            sm_metadata_fields[field_name] = []
            for ln in range(extent):
                field_line = check_call(
                    self.simulator_model.get_simulator_field_line, field, ln
                )
                sm_metadata_fields[field_name].append(field_line)

        return sm_metadata_fields

    @property
    def supported_units(self):
        try:
            supported_units = self.metadata["units"][0]
        except (KeyError, IndexError):
            raise KIMCalculatorError(
                "ERROR: Unable to determine supported units of "
                "simulator model {}.".format(self.model_name)
            )

        return supported_units

    @property
    def atom_style(self):
        """
        See if a 'model-init' field exists in the SM metadata and, if
        so, whether it contains any entries including an "atom_style"
        command.  This is specific to LAMMPS SMs and is only required
        for using the LAMMPSrun calculator because it uses
        lammps.inputwriter to create a data file.  All other content in
        'model-init', if it exists, is ignored.
        """
        atom_style = None
        for ln in self.metadata.get("model-init", []):
            if ln.find("atom_style") != -1:
                atom_style = ln.split()[1]

        return atom_style

    @property
    def model_defn(self):
        return self.metadata["model-defn"]

    def destroy(self):
        if self.initialized:
            kimpy.simulator_model.destroy(self.simulator_model)
            self.initialized = False
