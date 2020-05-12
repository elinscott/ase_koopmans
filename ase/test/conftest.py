import os
from pathlib import Path
from subprocess import Popen, PIPE

import pytest

from ase.utils import workdir
from ase.test.factories import (Factories, CalculatorInputs,
                                make_factory_fixture, get_testing_executables)


@pytest.fixture(scope='session')
def enabled_calculators():
    from ase.calculators.calculator import names as calculator_names
    all_names = set(calculator_names)
    names = set(os.environ.get('ASE_TEST_CALCULATORS', '').split())
    for name in names:
        assert name in all_names
    return names


@pytest.fixture(scope='session', autouse=True)
def monkeypatch_disabled_calculators(request, enabled_calculators):
    from ase.test.testsuite import disable_calculators, test_calculator_names
    from ase.calculators.calculator import names as calculator_names
    test_calculator_names += enabled_calculators
    disable_calculators([name for name in calculator_names
                         if name not in enabled_calculators])


# Backport of tmp_path fixture from pytest 3.9.
# We want to be compatible with pytest 3.3.2 and pytest-xdist 1.22.1.
# These are provided with Ubuntu 18.04.
# Current Debian stable uses a newer libraries, so that should be OK.
@pytest.fixture
def tmp_path(tmpdir):
    # Avoid trouble since tmpdir can be a py._path.local.LocalPath
    return Path(str(tmpdir))


@pytest.fixture(autouse=True)
def use_tmp_workdir(tmp_path):
    # Pytest can on some systems provide a Path from pathlib2.  Normalize:
    path = Path(str(tmp_path))
    with workdir(path, mkdir=True):
        yield tmp_path


@pytest.fixture(scope='session')
def tkinter():
    import tkinter
    try:
        tkinter.Tk()
    except tkinter.TclError as err:
        pytest.skip('no tkinter: {}'.format(err))


@pytest.fixture(scope='session')
def plt(tkinter):
    matplotlib = pytest.importorskip('matplotlib')
    matplotlib.use('Agg')

    import matplotlib.pyplot as plt
    return plt


@pytest.fixture
def figure(plt):
    fig = plt.figure()
    yield fig
    plt.close(fig)


@pytest.fixture(scope='session')
def psycopg2():
    return pytest.importorskip('psycopg2')


@pytest.fixture(scope='session')
def datafiles():
    try:
        import asetest
    except ImportError:
        return {}
    else:
        return asetest.datafiles.paths


@pytest.fixture(scope='session')
def configured_executables():
    return get_testing_executables()


@pytest.fixture(scope='session')
def factories(configured_executables, datafiles, enabled_calculators):
    return Factories(configured_executables, datafiles)


abinit_factory = make_factory_fixture('abinit')
cp2k_factory = make_factory_fixture('cp2k')
dftb_factory = make_factory_fixture('dftb')
espresso_factory = make_factory_fixture('espresso')
gpaw_factory = make_factory_fixture('gpaw')
octopus_factory = make_factory_fixture('octopus')
siesta_factory = make_factory_fixture('siesta')


@pytest.fixture
def factory(request, factories):
    name, kwargs = request.param
    factory = factories[name]
    return CalculatorInputs(factory, kwargs)


def pytest_generate_tests(metafunc):
    from ase.test.factories import parametrize_calculator_tests
    parametrize_calculator_tests(metafunc)

    if 'seed' in metafunc.fixturenames:
        seeds = metafunc.config.getoption('seed')
        if len(seeds) == 0:
            seeds = [0, 1]
        else:
            seeds = list(map(int, seeds))
        metafunc.parametrize('seed', seeds)


class CLI:
    def ase(self, args):
        if isinstance(args, str):
            import shlex
            args = shlex.split(args)

        proc = Popen(['ase', '-T'] + args,
                     stdout=PIPE, stdin=PIPE)
        stdout, _ = proc.communicate(b'')
        status = proc.wait()
        assert status == 0
        return stdout.decode('utf-8')

    def shell(self, command, calculator_name=None):
        from ase.test.testsuite import runshellcommand
        runshellcommand(command, calculator_name=calculator_name)


@pytest.fixture(scope='session')
def datadir():
    from ase.test.testsuite import datadir
    return datadir


@pytest.fixture(scope='session')
def asap3():
    asap3 = pytest.importorskip('asap3')
    return asap3


@pytest.fixture(scope='session')
def cli():
    return CLI()


@pytest.fixture(autouse=True)
def arbitrarily_seed_rng(request):
    # We want tests to not use global stuff such as np.random.seed().
    # But they do.
    #
    # So in lieu of (yet) fixing it, we reseed and unseed the random
    # state for every test.  That makes each test deterministic if it
    # uses random numbers without seeding, but also repairs the damage
    # done to global state if it did seed.
    #
    # In order not to generate all the same random numbers in every test,
    # we seed according to a kind of hash:
    import numpy as np
    module_name = request.module
    function_name = request.function.__name__
    seed = hash((module_name, function_name)) % 123456789
    # (We should really use the full qualified name of the test method.)
    state = np.random.get_state()
    np.random.seed(seed)
    yield
    np.random.set_state(state)


def pytest_addoption(parser):
    parser.addoption('--seed', action='append', default=[],
                     help='Add a seed for tests where random number generators'
                          ' are involved. This option can be applied more'
                          ' than once.')
