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
from ase.calculators.lammpslib import LAMMPSlib
from ase.calculators.lammpsrun import LAMMPS
from ase.calculators.lammps import convert
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
    # of its simulator and its units from its metadata
    this_is_a_KIM_MO = _is_portable_model(extended_kim_id)

    # If this is a KIM Portable Model (supports KIM API), return support through
    # a KIM-compliant simulator
    if this_is_a_KIM_MO:

        if simulator is None:   # Default
            simulator = 'kimmodel'

        supported_species = _get_portable_model_supported_species(extended_kim_id,
                simulator)

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
            parameters['masses'] = []
            for i, species in enumerate(supported_species):
                if species not in atomic_numbers:
                    raise KIMCalculatorError(
                        'Unknown element species {0}.'.format(species))
                massstr = str(atomic_masses[atomic_numbers[species]])
                parameters['masses'].append(str(i + 1) + " " + massstr)

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

    #######################################################
    # If we get to here, the model is a KIM Simulator Model
    #######################################################
    simulator_name, supported_species, supported_units, atom_style = _get_simulator_model_info(
            extended_kim_id
    )

    # determine simulator
    if simulator is None:
        if simulator_name == 'ASAP':
            simulator = 'asap'
        elif simulator_name == 'LAMMPS':
            simulator = 'lammpslib'

    if simulator_name == "ASAP":
        # Initialize KIM SM object
        from kimpy import simulator_models as kimsm
        ksm = kimsm.ksm_object(extended_kim_id=extended_kim_id)
        param_filenames = ksm.get_model_param_filenames()

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

    elif simulator_name == "LAMMPS":

        if simulator == 'lammpsrun':
            # check options
            msg = _check_conflict_options(
                options, lammpsrun_not_allowed_options, simulator)
            if msg is not None:
                raise KIMCalculatorError(msg)

            # Set up kim_init and kim_interactions lines
            parameters = _get_params_for_LAMMPS_calculator(extended_kim_id,
                    supported_units, supported_species, atom_style)

            # Return LAMMPS calculator
            return LAMMPS(**parameters,
                          specorder=supported_species, keep_tmp_files=debug)

        elif simulator == 'lammpslib':
            # check options
            msg = _check_conflict_options(
                options, lammpslib_not_allowed_options, simulator)
            if msg is not None:
                raise KIMCalculatorError(msg)

            # Set up LAMMPS header commands lookup table

            # This units command actually has no effect, but is necessary because
            # LAMMPSlib looks in the header lines for units in order to set them
            # internally
            model_init = ['units ' + supported_units + "\n"]

            model_init.append("kim_init {} {}\n".format(extended_kim_id,
                supported_units))
            model_init.append('atom_modify map array sort 0 0\n')

            # Assign atom types to species
            atom_types = {}
            for i_s, s in enumerate(supported_species):
                atom_types[s] = i_s + 1

            kim_interactions = ["kim_interactions {}".format((' ').join(supported_species))]

            # Return LAMMPSlib calculator
            return LAMMPSlib(lammps_header=model_init,
                             lammps_name=None,
                             lmpcmds=kim_interactions,
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


def _get_kim_api_libexec_path():
    # Find the location of the libexec path for the KIM API
    try:
        libexec_path = subprocess.check_output(
            ["pkg-config", "--variable=libexecdir", "libkim-api"],
            universal_newlines=True).strip().rstrip("/")
    except subprocess.CalledProcessError:
        raise KIMCalculatorError(
            'ERROR: Unable to obtain libexec-path for KIM API from pkg-config.')

    return libexec_path


def _is_portable_model(extended_kim_id):
    """
    Returns True if the model specified is a KIM Portable Model (if it is not, then it
    must be a KIM Simulator Model -- there are no other types of models in KIM)
    """
    # Define location of `kim-api-collections-management-info' utility
    libexec_path = _get_kim_api_libexec_path()
    kim_api_cm_info_util = os.path.join(libexec_path, "kim-api",
            "kim-api-collections-info")

    try:
        item_type = subprocess.check_output([kim_api_cm_info_util, "type", extended_kim_id],
                universal_newlines=True)
        item_type = item_type.rstrip()
    except subprocess.CalledProcessError:
        raise KIMCalculatorError(
                'ERROR: Unable to call kim-api-collections-info util to '
                'determine whether item is Portable Model or Simulator Model.')

    return (item_type == 'portableModel')


def _get_portable_model_supported_species(extended_kim_id, requested_simulator):
    '''
    Determines what species a KIM Portable Model or KIM Simulator Model supports
    '''
    if requested_simulator in ["kimmodel", "lammpsrun"]:
        calc = KIMModelCalculator(extended_kim_id)
        supported_species = list(calc.get_kim_model_supported_species())
        calc.__del__()

    elif requested_simulator == "asap":
        from asap3 import OpenKIMcalculator
        calc = OpenKIMcalculator(extended_kim_id)
        supported_species = list(calc.get_supported_elements())
        if hasattr(calc, '__del__'):
            calc.__del__()

    return supported_species


def _get_simulator_model_info(extended_kim_id):
    """
    Retrieve Simulator Model metadata including its native simulator, supported species,
    and units
    """
    # FIXME: This should be updated to use the appropriate bindings in kimpy
    libexec_path = _get_kim_api_libexec_path()
    kim_api_sm_util = os.path.join(libexec_path, "kim-api", "kim-api-simulator-model")

    try:
        sm_metadata = subprocess.check_output([kim_api_sm_util, extended_kim_id,
            "smspec-file", "data"], universal_newlines=True)
    except subprocess.CalledProcessError:
        raise KIMCalculatorError(
                'ERROR: Unable to call kim-api-collections-info util to '
                'retrieve Simulator Model metadata.')

    # Parse metadata for simulator-name
    simulator_name = re.search(r"\"simulator-name\"\s+\"([A-Za-z0-9]+)\"", sm_metadata)
    if simulator_name is None:
        raise KIMCalculatorError("ERROR: Unable to determine simulator name of "
                "item {}.".format(extended_kim_id))
    else:
        simulator_name = simulator_name.group(1)

    # Parse metadata for species
    supported_species = re.search(r"\"supported-species\"\s+\"([A-Za-z0-9\s]+)\"", sm_metadata)
    if supported_species is None:
        raise KIMCalculatorError("ERROR: Unable to determine supported species of "
                "item {}.".format(extended_kim_id))
    else:
        supported_species = supported_species.group(1).split()

    # Parse metadata for units
    supported_units = re.search(r"\"units\"\s+\"([A-Za-z0-9\s]+)\"", sm_metadata)
    if supported_units is None:
        raise KIMCalculatorError("ERROR: Unable to determine supported units of "
                "item {}.".format(extended_kim_id))
    else:
        supported_units = supported_units.group(1)

    return simulator_name, supported_species, supported_units


def _get_params_for_LAMMPS_calculator(extended_kim_id, supported_units,
        supported_species, atom_style):
    '''
    Extract parameters for LAMMPS calculator from model definition lines.
    Returns a dictionary with entries for "pair_style" and "pair_coeff".
    Expects there to be only one "pair_style" line. There can be multiple
    "pair_coeff" lines (result is returned as a list).
    '''
    parameters = {}

    # In case the SM supplied its own atom_style in its model-init -- only needed
    # because lammpsrun writes data files and needs to know the proper format
    parameters['atom_style'] = atom_style

    # Set units to prevent them from defaulting to metal
    parameters['units'] = supported_units

    parameters['model_init'] = ["kim_init {} {}\n".format(extended_kim_id,
        supported_units)]

    parameters['kim_interactions'] = "kim_interactions {}\n".format((' ').join(supported_species))

    # For every species in "supported_species", add an entry to the
    # "masses" key in dictionary "parameters".
    parameters['masses'] = []
    for i, species in enumerate(supported_species):
        if species not in atomic_numbers:
            raise KIMCalculatorError(
                'Unknown element species {}.'.format(species))
        massstr = str(convert(atomic_masses[atomic_numbers[species]], "mass", "ASE",
            supported_units))
        parameters['masses'].append(str(i + 1) + " " + massstr)

    return parameters


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
