"""
Knowledgebase of Interatomic Models (KIM) Calculator for ASE written by:

Ellad B. Tadmor
Mingjian Wen
University of Minnesota

This calculator selects an appropriate calculator for a KIM model depending on
whether it supports the KIM application programming interface (API) or is a
KIM Simulator Model. For more information on KIM, visit https://openkim.org.
"""

from __future__ import absolute_import
from __future__ import print_function
from __future__ import division
import re
import os
import subprocess
from ase.data import atomic_masses, atomic_numbers
try:
    from kimpy import simulator_models as kimsm
    kimsm_loaded = True
except Exception:
    kimsm_loaded = False
from ase.calculators.lammpslib import LAMMPSlib
from ase.calculators.lammpsrun import LAMMPS
from .kimmodel import KIMModelCalculator
from .exceptions import KIMCalculatorError


def KIM(extended_kim_id, simulator=None, options=None, debug=False):
    """Calculator for interatomic models archived in the Open Knowledgebase
    of Interatomic Models (OpenKIM) at https://openkim.org

    Parameters
    ----------
    extended_kim_id: string
        Extended KIM ID of the KIM interatomic model (for details see:
        https://openkim.org/about-kim-ids/)

    simulator: string (optional)
        Name of the simulator to be used. This is the name of the ASE calculator
        that will be used to run the KIM interatomic model. Supported simulators
        include ``kimmodel``, ``lammpslib``, ``lammpsrun`` and ``asap3``.
        If ``None``, simulator is determined automatically based on
        ``extended_kim_id``.

    options: dictionary (optional)
        Additional options passed to the initializer of the selected simulator.
        If the ``simulator="kimmodel"``, possible options are:

        options = {'neigh_skin_ratio': 0.2, 'release_GIL': False}

        where ``neigh_skin_ratio`` provides the skin (in percentage of cutoff)
        used to determine the neighbor list, and ``release_GIL`` determines
        whether to release the python GIL, which allows a KIM model to be run
        with multiple threads.

        See the LAMMPS calculators doc page
        https://wiki.fysik.dtu.dk/ase/ase/calculators/lammps.html
        for available options for ``lammpsrun`` and ``lammpslib``.

    debug: bool (optional)
        If ``True``, turn on the debug mode to output extra information.

    Note
    ----
    This calculator is actually a wrapper that returns a calculator based on the
    arguments ``extended_kim_id`` and ``simulator``.
    """

    # options set internally in this calculator
    kimmodel_not_allowed_options = ['modelname', 'debug']
    lammpsrun_not_allowed_options = ['parameters', 'files', 'specorder',
                                     'keep_tmp_files']
    lammpslib_not_allowed_options = ['lammps_header', 'lmpcmds',
                                     'atom_types', 'log_file', 'keep_alive']
    asap_kimmo_not_allowed_options = ['name', 'verbose']
    asap_simmo_not_allowed_options = ['Params']
    if options is None:
        options = dict()

    # Determine whether this is a standard KIM Portable Model or a KIM Simulator Model
    # and retrieve its supported species.  If it is a Simulator Model, also get the name
    # of its simulator
    this_is_a_KIM_MO, supported_species, simulator_name = _get_model_info(extended_kim_id, simulator)

    # If this is a KIM Portable Model (supports KIM API), return support through
    # a KIM-compliant simulator
    if this_is_a_KIM_MO:

        if simulator is None:   # Default
            simulator = 'kimmodel'

        if simulator == 'kimmodel':
            msg = _check_conflict_options(
                options, kimmodel_not_allowed_options, simulator)
            if msg is not None:
                raise KIMCalculatorError(msg)
            else:
                return KIMModelCalculator(extended_kim_id, debug=debug,
                                          **options)

        elif simulator == 'asap':
            try:
                from asap3 import OpenKIMcalculator
            except ImportError as e:
                raise ImportError(str(e) + ' You need to install asap3 first.')

            msg = _check_conflict_options(
                options, asap_kimmo_not_allowed_options, simulator)
            if msg is not None:
                raise KIMCalculatorError(msg)
            else:
                return OpenKIMcalculator(name=extended_kim_id, verbose=debug,
                                         **options)

        elif simulator == 'lammpsrun':
            msg = _check_conflict_options(
                options, lammpsrun_not_allowed_options, simulator)
            if msg is not None:
                raise KIMCalculatorError(msg)

            param_filenames = []  # no parameter files to pass
            parameters = {}
            parameters['pair_style'] = 'kim ' + \
                extended_kim_id.strip() + os.linesep
            parameters['pair_coeff'] = [
                '* * ' + ' '.join(supported_species) + os.linesep]
            parameters['model_init'] = []
            parameters['model_post'] = []
            parameters['mass'] = []
            for i, species in enumerate(supported_species):
                if species not in atomic_numbers:
                    raise KIMCalculatorError(
                        'Unknown element species {0}.'.format(species))
                massstr = str(atomic_masses[atomic_numbers[species]])
                parameters['mass'].append(str(i + 1) + " " + massstr)

            # Return LAMMPS calculator
            return LAMMPS(**parameters, files=param_filenames,
                          specorder=supported_species, keep_tmp_files=debug)

        elif simulator == 'lammpslib':
            raise KIMCalculatorError(
                '"lammpslib" does not support KIM model. try "lammpsrun".')
        else:
            raise KIMCalculatorError(
                'Unsupported simulator "{}" requested to run KIM Models.'
                .format(simulator))

    # If we get to here, the model is a KIM Simulator Model

    # Initialize KIM SM object
    ksm = kimsm.ksm_object(extended_kim_id=extended_kim_id)
    param_filenames = ksm.get_model_param_filenames()

    # Double check that the extended KIM ID of the Simulator Model
    # matches the expected value. (If not, the KIM SM is corrupted.)
    SM_extended_kim_id = ksm.get_model_extended_kim_id()
    if extended_kim_id != SM_extended_kim_id:
        raise KIMCalculatorError(
            'SM extended KIM ID ("{}") does not match expected value '
            ' ("{}").'.format(SM_extended_kim_id, extended_kim_id))

    # determine simulator
    if simulator is None:
        if simulator_name == 'asap':
            simulator = 'asap'
        elif simulator_name == 'lammps':
            simulator = 'lammpslib'

    #  Get model definition from SM metadata
    model_defn = ksm.get_model_defn_lines()
    if len(model_defn) == 0:
        raise KIMCalculatorError(
            'model-defn is an empty list in metadata file of '
            'Simulator Model "{}".'.format(extended_kim_id))
    if "" in model_defn:
        raise KIMCalculatorError(
            'model-defn contains one or more empty strings in metadata '
            'file of Simulator Model "{}".'.format(extended_kim_id))

    if simulator_name == "asap":
        try:
            from asap3 import EMT, EMTMetalGlassParameters, EMTRasmussenParameters
        except ImportError as e:
            raise ImportError(str(e) + ' You need to install asap3 first.')

        # check options
        msg = _check_conflict_options(
            options, asap_simmo_not_allowed_options, simulator)
        if msg is not None:
            raise KIMCalculatorError(msg)

        # Verify units (ASAP models are expected to work with "ase" units)
        supported_units = ksm.get_model_units().lower().strip()
        if supported_units != "ase":
            raise KIMCalculatorError(
                'KIM Simulator Model units are "{}", but expected to '
                'be "ase" for ASAP.'.format(supported_units))

        # There should be only one model_defn line
        if len(model_defn) != 1:
            raise KIMCalculatorError(
                'model-defn contains {} lines, but should only contain '
                'one line for an ASAP model.'.format(len(model_defn)))

        # Return calculator
        unknown_potential = False
        if model_defn[0].lower().strip().startswith("emt"):
            # pull out potential parameters
            pp = ''
            mobj = re.search(r"\(([A-Za-z0-9_\(\)]+)\)", model_defn[0])
            if mobj is not None:
                pp = mobj.group(1).strip().lower()
            if pp == '':
                calc = EMT()
            elif pp.startswith('emtrasmussenparameters'):
                calc = EMT(Params=EMTRasmussenParameters())
            elif pp.startswith('emtmetalglassparameters'):
                calc = EMT(Params=EMTMetalGlassParameters())
            else:
                unknown_potential = True

        if unknown_potential:
            raise KIMCalculatorError(
                'Unknown model "{}" for simulator ASAP.'.format(model_defn[0]))
        else:
            calc.set_subtractE0(False)  # Use undocumented feature for the EMT
            # calculators to take the energy of an
            # isolated atoms as zero. (Otherwise it
            # is taken to be that of perfect FCC.)
            return calc

    elif simulator_name == "lammps":

        param_filenames_for_lammps = list(param_filenames)
        if simulator == 'lammpsrun':
            # Remove path from parameter file names since lammpsrun copies all
            # files into a tmp directory, so path should not appear on
            # in LAMMPS commands
            param_filenames_for_lammps = [os.path.basename(i)
                                          for i in param_filenames_for_lammps]

        # Build atom species and type lists based on all supported species.
        # This means that the LAMMPS simulation will be defined to have
        # as many atom types as are supported by the SM and each atom will
        # be assigned a type based on its species (in the order that the
        # species are defined in the SM).
        supported_species = ksm.get_model_supported_species()
        atom_type_sym_list_string = ' '.join(supported_species)
        atom_type_num_list_string = ' '.join(
            [str(atomic_numbers[s]) for s in supported_species])

        # Process KIM templates in model_defn lines
        for i in range(0, len(model_defn)):
            model_defn[i] = kimsm.template_substitution(
                model_defn[i], param_filenames_for_lammps, ksm.sm_dirname,
                atom_type_sym_list_string, atom_type_num_list_string)

        # Get model init lines
        model_init = ksm.get_model_init_lines()

        # Process KIM templates in model_init lines
        for i in range(0, len(model_init)):
            model_init[i] = kimsm.template_substitution(
                model_init[i], param_filenames_for_lammps, ksm.sm_dirname,
                atom_type_sym_list_string, atom_type_num_list_string)

        # Get model supported units
        supported_units = ksm.get_model_units().lower().strip()

        if simulator == 'lammpsrun':
            # check options
            msg = _check_conflict_options(
                options, lammpsrun_not_allowed_options, simulator)
            if msg is not None:
                raise KIMCalculatorError(msg)

            # add cross-platform line separation to model definition lines
            model_defn = [s + os.linesep for s in model_defn]

            # Extract parameters for calculator from model definition lines
            parameters = _get_params_for_LAMMPS_calculator(model_defn,
                                                           supported_species)

            # Add units to parameters
            parameters["units"] = supported_units

            # add cross-platform line separation to model definition lines
            model_init = [s + os.linesep for s in model_init]

            # Add init lines to parameter list
            _add_init_lines_to_parameters(parameters, model_init)

            # Return LAMMPS calculator
            return LAMMPS(**parameters, files=param_filenames,
                          specorder=supported_species, keep_tmp_files=debug)

        elif simulator == 'lammpslib':
            # check options
            msg = _check_conflict_options(
                options, lammpslib_not_allowed_options, simulator)
            if msg is not None:
                raise KIMCalculatorError(msg)

            # Set up LAMMPS header commands lookup table
            model_init.insert(0, 'atom_modify map array sort 0 0')
            if not any("atom_style" in s.lower() for s in model_init):
                model_init.insert(0, 'atom_style atomic')
            model_init.insert(
                0, 'units ' + supported_units.strip())     # units

            # Assign atom types to species
            atom_types = {}
            for i_s, s in enumerate(supported_species):
                atom_types[s] = i_s + 1

            # Return LAMMPSlib calculator
            return LAMMPSlib(lammps_header=model_init,
                             lammps_name=None,
                             lmpcmds=model_defn,
                             atom_types=atom_types,
                             log_file='lammps.log',
                             keep_alive=True,
                             **options)

        else:
            raise KIMCalculatorError(
                'Unknown LAMMPS calculator: "{}".'.format(simulator))

    else:
        raise KIMCalculatorError(
            'Unsupported simulator: "{}".'.format(simulator_name))


def _get_model_info(extended_kim_id, requested_simulator):
    '''
    Determine whether a model corresponds to either a Portable Model or Simulator Model
    and what species it supports.  If it is a Simulator Model, also return the name of
    its simulator.
    '''
    # Find the location of the `kim-api-collections-management-info' utility
    try:
        libexec_path = subprocess.check_output(
            ["pkg-config", "--variable=libexecdir", "libkim-api"],
            universal_newlines=True).strip().rstrip("/")
    except subprocess.CalledProcessError:
        raise KIMCalculatorError(
            'ERROR: Unable to obtain libexec-path for KIM API from pkg-config.')

    kim_api_cm_info_util = os.path.join(libexec_path, "kim-api",
            "kim-api-collections-info")

    # Determine whether this is a Portable Model or Simulator Model
    try:
        item_type = subprocess.check_output([kim_api_cm_info_util, "type", extended_kim_id],
                universal_newlines=True)
        item_type = item_type.rstrip()
    except subprocess.CalledProcessError:
        raise KIMCalculatorError(
                'ERROR: Unable to call kim-api-collections-info util to '
                'determine whether item is Portable Model or Simulator Model.')

    if item_type == 'portableModel':
        this_is_a_KIM_MO = True
        simulator_name = None

        if requested_simulator == "kimmodel":
            with KIMModelCalculator(extended_kim_id) as calc:
                supported_species = list(calc.get_kim_model_supported_species())

        elif requested_simulator == "asap":
            from asap3 import OpenKIMcalculator
            with OpenKIMcalculator(extended_kim_id) as calc:
                supported_species = list(calc.get_supported_elements())

    elif item_type == 'simulatorModel':
        this_is_a_KIM_MO = False

        # Retrieve Simulator Model metadata
        try:
            sm_metadata = subprocess.check_output([kim_api_cm_info_util, extended_kim_id,
                "smspec-file", "data"], universal_newlines=True)
        except subprocess.CalledProcessError:
            raise KIMCalculatorError(
                    'ERROR: Unable to call kim-api-collections-info util to '
                    'retrieve Simulator Model metadata.')

        # Parse metadata for simulator-name
        simulator_name = re.search("\"simulator-name\"\s+\"([A-Za-z0-9]+)\"", sm_metadata)
        if simulator_name is None:
            raise KIMCalculatorError("ERROR: Unable to determine simulator name of "
                    "item {}.".format(extended_kim_id))
        else:
            simulator_name = simulator_name.groups(1)

        # Parse metadata for species
        supported_species = re.search("\"supported-species\"\s+\"([A-Za-z0-9\s]+)\"", sm_metadata)
        if supported_species is None:
            raise KIMCalculatorError("ERROR: Unable to determine supported species of "
                    "item {}.".format(extended_kim_id))
        else:
            supported_species = supported_species.groups(1)
    else:
        raise KIMCalculatorError("ERROR: Item {} has type {} and is not a Portable "
            "Model or Simulator Model.".format(extended_kim_id, item_type))

    return this_is_a_KIM_MO, supported_species, simulator_name


def _get_params_for_LAMMPS_calculator(model_defn, supported_species):
    '''
    Extract parameters for LAMMPS calculator from model definition lines.
    Returns a dictionary with entries for "pair_style" and "pair_coeff".
    Expects there to be only one "pair_style" line. There can be multiple
    "pair_coeff" lines (result is returned as a list).
    '''
    parameters = {}
    parameters['pair_style'] = ''
    parameters['pair_coeff'] = []
    parameters['model_post'] = []
    found_pair_style = False
    found_pair_coeff = False
    for i in range(0, len(model_defn)):
        c = model_defn[i]
        if c.lower().startswith('pair_style'):
            if found_pair_style:
                raise KIMCalculatorError(
                    'More than one pair_style in metadata file.')
            found_pair_style = True
            parameters['pair_style'] = c.split(" ", 1)[1]
        elif c.lower().startswith('pair_coeff'):
            found_pair_coeff = True
            parameters['pair_coeff'].append(c.split(" ", 1)[1])
        else:
            parameters['model_post'].append(c)
    if not found_pair_style:
        raise KIMCalculatorError('pair_style not found in metadata file.')
    if not found_pair_coeff:
        raise KIMCalculatorError('pair_coeff not found in metadata file.')

    # For every species in "supported_species", add an entry to the
    # "mass" key in dictionary "parameters".
    parameters['mass'] = []
    for i, species in enumerate(supported_species):
        if species not in atomic_numbers:
            raise KIMCalculatorError(
                'Unknown element species {0}.'.format(species))
        massstr = str(atomic_masses[atomic_numbers[species]])
        parameters['mass'].append(str(i + 1) + " " + massstr)

    return parameters


def _add_init_lines_to_parameters(parameters, model_init):
    '''
    Add Simulator Model initialization lines to the parameter list for LAMMPS
    if there are any.
    '''
    parameters['model_init'] = []
    for i in range(0, len(model_init)):
        parameters['model_init'].append(model_init[i])


def _check_conflict_options(options, not_allowed_options, simulator):
    """Check whether options is in not_allowed options"""
    s1 = set(options)
    s2 = set(not_allowed_options)
    common = s1.intersection(s2)
    if common:
        msg1 = 'Simulator "{}" does not support argument(s): '.format(
            simulator)
        msg2 = ', '.join(['"{}"'.format(s) for s in common])
        msg3 = ' provided in "options", because it is (they are) determined '
        msg4 = 'internally within the KIM calculator.'
        return msg1 + msg2 + msg3 + msg4
    else:
        msg = None
    return msg
