"""
Knowledgebase of Interatomic Models (KIM) Calculator for ASE written by:

Ellad B. Tadmor
Mingjian Wen
Daniel S. Karls
University of Minnesota

This calculator functions as a wrapper that selects an appropriate calculator for a
given KIM model depending on whether it supports the KIM application programming
interface (API) or not. For more information on KIM, visit https://openkim.org.
"""

from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

import re
import os

from ase.data import atomic_masses, atomic_numbers
from ase.calculators.lammpslib import LAMMPSlib
from ase.calculators.lammpsrun import LAMMPS
from ase.calculators.lammps import convert

from .kimmodel import KIMModelCalculator

try:
    import kimpy
except ImportError:
    raise RuntimeError("kimpy not found; KIM calculator will not work")


class KIMCalculatorError(Exception):
    pass


def KIM(extended_kim_id, simulator=None, options=None, debug=False):
    """Calculator wrapper for OpenKIM models

    Returns a suitable calculator that can be used with any model
    archived in the Open Knowledgebase of Interatomic Models (OpenKIM)
    at https://openkim.org.  There are two kinds of models in KIM:
    Portable Models (PMs), which contain an '__MO_' in their name and
    can be used with any KIM-compliant simulator, and Simulator Models
    (SMs), which have an '__SM_' in their name and are essentially just
    wrappers around native commands in a specific simulator.  In
    general, 'simulator' can be omitted, in which case this function
    will automatically determine a calculator compatible with the
    specified model and return it.

    Parameters
    ----------
    extended_kim_id : str
        Extended KIM ID of the interatomic model (for details, see
        https://openkim.org/doc/schema/kim-ids)

    simulator : str, optional
        Name of the ASE calculator that will be used.  Supported
        simulators currently include 'kimmodel', 'lammpslib',
        'lammpsrun' and 'asap'.  However, not all of these values are
        compatible with all KIM models; for example, a KIM SM created
        for LAMMPS must either use either 'lammpslib' or 'lammpsrun'.
        If None, a compatible simulator is determined automatically
        based on `extended_kim_id`.

    options : dict, optional
        Additional options passed to the initializer of the selected
        simulator.  If `simulator`=='kimmodel', possible options are:

        options = {'neigh_skin_ratio': 0.2, 'release_GIL': False}

        where 'neigh_skin_ratio' provides the skin (as a factor of the
        model cutoff) used to determine the neighbor list, and
        'release_GIL' determines whether to release the python GIL,
        which allows a KIM model to be run with multiple threads.

        See the ASE LAMMPS calculators doc page
        (https://wiki.fysik.dtu.dk/ase/ase/calculators/lammps.html) for
        available options for the lammpslib and lammpsrun calculators.

    debug: bool, optional
        If True, detailed information is printed to stdout.  If the
        lammpsrun calculator is being used, this also serves as the
        value of the `keep_tmp_files` option.

    Returns
    -------
    calc : Calculator
        An ASE calculator.  Currently, this will be an instance of
        KIMModelCalculator, LAMMPS (the lammpsrun calculator), or
        LAMMPSlib, which are all defined in the ASE codebase, or an
        instance of one of either OpenKIMcalculator or EMT, which are
        defined in the asap3 codebase.

    Raises
    ------
    KIMCalculatorError
        Blanket exception type used to handle errors that arise related
        to using incompatible combinations of values for the arguments
        or from errors produced by kimpy

    See Also
    --------
        asap3.Internal.OpenKIMcalculator.OpenKIMcalculator
        asap3.Internal.BuiltinPotentials.EMT
        ase.calculators.kim.kimmodel.KIMModelCalculator
        ase.calculators.lammpslib.LAMMPSlib
        ase.calculators.lammpsrun.LAMMPS
    """

    # calculator to return
    calc = None

    # options set internally in this calculator
    kimmodel_not_allowed_options = ["modelname", "debug"]
    lammpsrun_not_allowed_options = [
        "parameters",
        "files",
        "specorder",
        "keep_tmp_files",
    ]
    lammpslib_not_allowed_options = [
        "lammps_header",
        "lmpcmds",
        "atom_types",
        "log_file",
        "keep_alive",
    ]
    asap_kimpm_not_allowed_options = ["name", "verbose"]
    asap_kimsm_not_allowed_options = ["Params"]
    if options is None:
        options = dict()

    # If this is a KIM Portable Model (supports KIM API), return support through
    # a KIM-compliant simulator
    if _is_portable_model(extended_kim_id):
        if simulator is None:  # Default
            simulator = "kimmodel"

        if simulator == "kimmodel":
            _check_conflict_options(options, kimmodel_not_allowed_options, simulator)
            calc = KIMModelCalculator(extended_kim_id, debug=debug, **options)
            return calc

        elif simulator == "asap":
            try:
                from asap3 import OpenKIMcalculator
            except ImportError as e:
                raise ImportError(str(e) + " You need to install asap3 first.")

            _check_conflict_options(options, asap_kimpm_not_allowed_options, simulator)
            calc = OpenKIMcalculator(name=extended_kim_id, verbose=debug, **options)
            return calc

        elif simulator == "lammpsrun":

            _check_conflict_options(options, lammpsrun_not_allowed_options, simulator)

            supported_species = _get_kim_pm_supported_species(extended_kim_id)

            # Set up kim_init and kim_interactions lines
            parameters = _get_params_for_LAMMPS_calculator(
                extended_kim_id,
                supported_units="metal",
                supported_species=supported_species,
                atom_style=None,
            )

            # Return LAMMPS calculator
            calc = LAMMPS(
                **parameters,
                specorder=supported_species,
                keep_tmp_files=debug,
                **options
            )
            return calc

        elif simulator == "lammpslib":
            raise KIMCalculatorError(
                '"lammpslib" calculator does not support KIM Portable Models. Try '
                'using the "lammpsrun" calculator.'
            )
        else:
            raise KIMCalculatorError(
                'Unsupported simulator "{}" requested to run KIM Portable Model.'.format(
                    simulator
                )
            )

    #######################################################
    # If we get to here, the model is a KIM Simulator Model
    #######################################################
    (
        simulator_name,
        supported_species,
        supported_units,
        model_defn,
        atom_style,
    ) = _get_simulator_model_info(extended_kim_id)

    # Handle default behavior for 'simulator'
    if simulator is None:
        if simulator_name == "ASAP":
            simulator = "asap"
        elif simulator_name == "LAMMPS":
            simulator = "lammpslib"

    if simulator_name == "ASAP":
        # check options
        _check_conflict_options(options, asap_kimsm_not_allowed_options, simulator)

        return _asap_kimsm_calculator(extended_kim_id, model_defn, supported_units)

    elif simulator_name == "LAMMPS":

        if simulator == "lammpsrun":
            # check options
            _check_conflict_options(options, lammpsrun_not_allowed_options, simulator)

            # Set up kim_init and kim_interactions lines
            parameters = _get_params_for_LAMMPS_calculator(
                extended_kim_id, supported_units, supported_species, atom_style
            )

            # Return LAMMPS calculator
            calc = LAMMPS(
                **parameters, specorder=supported_species, keep_tmp_files=debug
            )
            return calc

        elif simulator == "lammpslib":
            # check options
            _check_conflict_options(options, lammpslib_not_allowed_options, simulator)

            # Set up LAMMPS header commands lookup table

            # This units command actually has no effect, but is necessary because
            # LAMMPSlib looks in the header lines for units in order to set them
            # internally
            model_init = ["units " + supported_units + os.linesep]

            model_init.append(
                "kim_init {} {}{}".format(extended_kim_id, supported_units, os.linesep)
            )
            model_init.append("atom_modify map array sort 0 0" + os.linesep)

            # Assign atom types to species
            atom_types = {}
            for i_s, s in enumerate(supported_species):
                atom_types[s] = i_s + 1

            kim_interactions = [
                "kim_interactions {}".format((" ").join(supported_species))
            ]

            # Return LAMMPSlib calculator
            calc = LAMMPSlib(
                lammps_header=model_init,
                lammps_name=None,
                lmpcmds=kim_interactions,
                atom_types=atom_types,
                log_file="lammps.log",
                keep_alive=True,
                **options
            )
            return calc

        else:
            raise KIMCalculatorError(
                'Unknown LAMMPS calculator: "{}".'.format(simulator)
            )

    else:
        raise KIMCalculatorError('Unsupported simulator: "{}".'.format(simulator_name))


def _is_portable_model(extended_kim_id):
    """
    Returns True if the model specified is a KIM Portable Model (if it
    is not, then it must be a KIM Simulator Model -- there are no other
    types of models in KIM)
    """
    col = check_call(kimpy.collections.create)

    model_type = check_call(col.get_item_type, extended_kim_id)

    kimpy.collections.destroy(col)

    return model_type == kimpy.collection_item_type.portableModel


def _get_simulator_model_info(extended_kim_id):
    """
    Retrieve Simulator Model metadata including its native simulator,
    supported species, and units
    """
    # Create a KIM API simulator Model object for this model
    kim_simulator_model = check_call(kimpy.simulator_model.create, extended_kim_id)

    # Retrieve simulator name (disregard simulator version)
    simulator_name, _ = kim_simulator_model.get_simulator_name_and_version()

    # Retrieve supported species
    num_supported_species = kim_simulator_model.get_number_of_supported_species()
    if num_supported_species == 0:
        raise KIMCalculatorError(
            "ERROR: Unable to determine supported species of "
            "simulator model {}.".format(extended_kim_id)
        )

    supported_species = []
    for spec_code in range(num_supported_species):
        species = check_call(kim_simulator_model.get_supported_species, spec_code)
        supported_species.append(species)

    # Need to close template map to access simulator model metadata
    kim_simulator_model.close_template_map()

    # Retrieve simulator model metadata
    sm_metadata_fields = {}
    num_metadata_fields = kim_simulator_model.get_number_of_simulator_fields()
    for field in range(num_metadata_fields):
        extent, field_name = check_call(
            kim_simulator_model.get_simulator_field_metadata, field
        )
        sm_metadata_fields[field_name] = []
        for ln in range(extent):
            field_line = check_call(
                kim_simulator_model.get_simulator_field_line, field, ln
            )
            sm_metadata_fields[field_name].append(field_line)

    # Grab units from simulator model metadata
    try:
        supported_units = sm_metadata_fields["units"][0]
    except (KeyError, IndexError):
        raise KIMCalculatorError(
            "ERROR: Unable to determine supported units of "
            "simulator model {}.".format(extended_kim_id)
        )

    # See if a 'model-init' field that contains an "atom_style" command is listed in
    # the simulator model metadata.  This is specific to LAMMPS SMs and is only required
    # for using the LAMMPSrun calculator because it uses lammps.inputwriter to create a
    # data file.  All other content in 'model-init', if it exists, is ignored
    atom_style = None
    try:
        for ln in sm_metadata_fields["model-init"]:
            if ln.find("atom_style") != -1:
                atom_style = ln.split()[1]
    except KeyError:
        pass

    # Clean up KIM API Simulator Model object
    kimpy.simulator_model.destroy(kim_simulator_model)

    return (
        simulator_name,
        tuple(supported_species),
        supported_units,
        sm_metadata_fields["model-defn"],
        atom_style,
    )


def _get_kim_pm_supported_species(extended_kim_id):
    """Gets species supported by a KIM Portable Model"""
    with KIMModelCalculator(extended_kim_id) as kim_calc:
        supported_species, _ = kim_calc.get_model_supported_species_and_codes()

    return tuple(supported_species)


def get_model_supported_species(extended_kim_id):
    """Convenience function for simulator codes"""
    if _is_portable_model(extended_kim_id):
        supported_species = _get_kim_pm_supported_species(extended_kim_id)
    else:
        _, supported_species, _, _, _ = _get_simulator_model_info(extended_kim_id)

    return supported_species


def _asap_kimsm_calculator(extended_kim_id, model_defn, supported_units):
    """
    Given the 'model-defn' entry from the metadata file (smspec.edn) of
    an ASAP KIM SM, construct an applicable asap3 calculator.  The
    supported units read from the metadata file are simply checked to
    ensure they correspond to ASE's native units.
    """
    # Check model_defn to make sure there's only one element in it that is a non-empty
    # string
    if len(model_defn) == 0:
        raise KIMCalculatorError(
            "model-defn is an empty list in metadata file of Simulator Model {}"
            "".format(extended_kim_id)
        )
    elif len(model_defn) > 1:
        raise KIMCalculatorError(
            "model-defn should contain only one entry for an ASAP model (found {} "
            "lines)".format(len(model_defn))
        )

    if "" in model_defn:
        raise KIMCalculatorError(
            "model-defn contains an empty string in metadata file of Simulator Model {}"
            "".format(extended_kim_id)
        )

    model_defn = model_defn[0].strip().lower()

    try:
        from asap3 import EMT, EMTMetalGlassParameters, EMTRasmussenParameters
    except ImportError as e:
        raise ImportError(str(e) + " You need to install asap3 first.")

    # Verify units (ASAP models are expected to work with "ase" units)
    if supported_units != "ase":
        raise KIMCalculatorError(
            'KIM Simulator Model units are "{}", but expected to '
            'be "ase" for ASAP.'.format(supported_units)
        )

    # Return calculator
    if model_defn.startswith("emt"):
        # pull out potential parameters
        mobj = re.search(r"\(([a-z0-9_\(\)]+)\)", model_defn)
        if mobj is None:
            asap_calc = EMT()
        else:
            pp = mobj.group(1)

            if pp.startswith("emtrasmussenparameters"):
                asap_calc = EMT(parameters=EMTRasmussenParameters())
            elif pp.startswith("emtmetalglassparameters"):
                asap_calc = EMT(parameters=EMTMetalGlassParameters())
            else:
                raise KIMCalculatorError(
                    'Unknown model "{}" for simulator ASAP.'.format(model_defn)
                )

    # Use undocumented feature for the EMT asap_calculators to take the energy of an
    # isolated atoms as zero. (Otherwise it is taken to be that of perfect FCC.)
    asap_calc.set_subtractE0(False)

    return asap_calc


def _get_params_for_LAMMPS_calculator(
    extended_kim_id, supported_units, supported_species, atom_style
):
    """
    Extract parameters for LAMMPS calculator from model definition lines.
    Returns a dictionary with entries for "pair_style" and "pair_coeff".
    Expects there to be only one "pair_style" line. There can be multiple
    "pair_coeff" lines (result is returned as a list).
    """
    parameters = {}

    # In case the SM supplied its own atom_style in its model-init -- only needed
    # because lammpsrun writes data files and needs to know the proper format
    if atom_style:
        parameters["atom_style"] = atom_style

    # Set units to prevent them from defaulting to metal
    parameters["units"] = supported_units

    parameters["model_init"] = [
        "kim_init {} {}{}".format(extended_kim_id, supported_units, os.linesep)
    ]

    parameters["kim_interactions"] = "kim_interactions {}{}".format(
        (" ").join(supported_species), os.linesep
    )

    # For every species in "supported_species", add an entry to the
    # "masses" key in dictionary "parameters".
    parameters["masses"] = []
    for i, species in enumerate(supported_species):
        if species not in atomic_numbers:
            raise KIMCalculatorError("Unknown element species {}.".format(species))
        massstr = str(
            convert(
                atomic_masses[atomic_numbers[species]], "mass", "ASE", supported_units
            )
        )
        parameters["masses"].append(str(i + 1) + " " + massstr)

    return parameters


def _check_error(error, msg):
    if error != 0 and error is not None:
        raise KIMCalculatorError('Calling "{}" failed.'.format(msg))


def check_call(f, *args):
    """
    Wrapper for functions that checks error codes, since many of the
    functions calls are actually python bindings to functions written in
    other languages.  Functions are assumed to return either an integer
    error code or a tuple whose last element is an integer error code.
    """
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


def _check_conflict_options(options, not_allowed_options, simulator):
    """Check whether options is in not_allowed options"""
    s1 = set(options)
    s2 = set(not_allowed_options)
    common = s1.intersection(s2)

    if common:
        options_in_not_allowed = ", ".join(['"{}"'.format(s) for s in common])

        msg = (
            'Simulator "{}" does not support argument(s): {} provided in "options", '
            "because it is (they are) determined internally within the KIM "
            "calculator".format(simulator, options_in_not_allowed)
        )

        raise KIMCalculatorError(msg)
