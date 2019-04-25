import os
import shutil
import importlib
from ase.calculators.calculator import names

builtins = {'eam', 'emt', 'lj', 'morse', 'tip3p'}

default_executables = {'abinit': ['abinit'],
                       'cp2k': ['cp2k_shell'],
                       'elk': ['elk', 'elk-lapw'],
                       'espresso': ['pw.x'],
                       'gromacs': ['gmx', 'gmx_d', 'gmx_mpi', 'gmx_mpi_d'],
                       'lammpsrun': ['lammps'],
                       'mopac': ['mopac', 'run_mopac7'], # run_mopac7: debian
                       'nwchem': ['nwchem'],
                       'octopus': ['octopus'],
                       'openmx': ['openmx'],
                       'siesta': ['siesta'],
                       }

python_modules = {'gpaw': 'gpaw',
                  'asap': 'asap3'}

def get_executable_env_var(name):
    return 'ASE_{}_COMMAND'.format(name.upper())


def detect(name):
    assert name in names
    d = {'name': name}

    if name in builtins:
        d['type'] = 'builtin'
        return d

    if name in python_modules:
        loader = importlib.find_loader(python_modules[name])
        if loader is not None:
            d['type'] = 'python'
            d['module'] = python_modules[name]
            d['path'] = loader.get_filename()
            return d

    envvar = get_executable_env_var(name)
    if envvar in os.environ:
        d['command'] = os.environ[envvar]
        d['envvar'] = envvar
        d['type'] = 'environment'
        return d

    if name in default_executables:
        commands = default_executables[name]
        for command in commands:
            fullpath = shutil.which(command)
            if fullpath:
                d['command'] = command
                d['fullpath'] = fullpath
                d['type'] = 'which'
                return d


def detect_calculators():
    configs = {}
    for name in names:
        result = detect(name)
        if result:
            configs[name] = result
    return configs


def format_configs(configs):
    messages = []
    for name in names:
        config = configs.get(name)

        if config is None:
            state = 'no'
        else:
            type = config['type']
            if type == 'builtin':
                state = 'yes, builtin: module ase.calculators.{name}'
            elif type == 'python':
                state = 'yes, python: {module} ▶ {path}'
            elif type == 'which':
                state = 'yes, shell command: {command} ▶ {fullpath}'
            else:
                state = 'yes, environment: ${envvar} ▶ {command}'

            state = state.format(**config)

        messages.append('{:<10s} {}'.format(name, state))
    return messages
