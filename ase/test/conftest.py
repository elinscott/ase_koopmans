import os
import pytest
from ase.utils import workdir
from pathlib import Path


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
    matplotlib.use('Agg', warn=False)

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
def executables():
    import json
    path = os.environ.get('ASE_EXECUTABLE_CONFIGFILE')
    if path is None:
        path = Path.home() / '.ase' / 'executables.json'
    path = Path(path)
    if path.exists():
        dct = json.loads(path.read_text())
    else:
        dct = {}
    return dct


factory_classes = {}


def factory(name):
    def deco(cls):
        cls.name = name
        factory_classes[name] = cls
        return cls
    return deco


@factory('abinit')
class AbinitFactory:
    def __init__(self, executable, pp_paths):
        self.executable = executable
        self.pp_paths = pp_paths

    def _base_kw(self):
        command = '{} < PREFIX.files > PREFIX.log'.format(self.executable)
        return dict(command=command,
                    pp_paths=self.pp_paths,
                    ecut=150,
                    chksymbreak=0,
                    toldfe=1e-3)

    def calc(self, **kwargs):
        from ase.calculators.abinit import Abinit
        kw = self._base_kw()
        kw.update(kwargs)
        return Abinit(**kw)

    @classmethod
    def fromconfig(cls, config):
        return AbinitFactory(config.executables['abinit'],
                             config.datafiles['abinit'])


@factory('cp2k')
class CP2KFactory:
    def __init__(self, executable):
        self.executable = executable

    def calc(self, **kwargs):
        from ase.calculators.cp2k import CP2K
        return CP2K(command=self.executable, **kwargs)

    @classmethod
    def fromconfig(cls, config):
        return CP2KFactory(config.executables['cp2k'])


@factory('espresso')
class EspressoFactory:
    def __init__(self, executable, pseudo_dir):
        self.executable = executable
        self.pseudo_dir = pseudo_dir

    def _base_kw(self):
        from ase.units import Ry
        return dict(ecutwfc=300 / Ry)

    def calc(self, **kwargs):
        from ase.calculators.espresso import Espresso
        command = '{} -in PREFIX.pwi > PREFIX.pwo'.format(self.executable)
        pseudopotentials = {}
        for path in self.pseudo_dir.glob('*.UPF'):
            fname = path.name
            # Names are e.g. si_lda_v1.uspp.F.UPF
            symbol = fname.split('_', 1)[0].capitalize()
            pseudopotentials[symbol] = fname

        kw = self._base_kw()
        kw.update(kwargs)
        return Espresso(command=command, pseudo_dir=str(self.pseudo_dir),
                        pseudopotentials=pseudopotentials,
                        **kw)

    @classmethod
    def fromconfig(cls, config):
        paths = config.datafiles['espresso']
        assert len(paths) == 1
        return cls(config.executables['espresso'], paths[0])


@factory('siesta')
class SiestaFactory:
    def __init__(self, executable, pseudo_path):
        self.executable = executable
        self.pseudo_path = pseudo_path

    def calc(self, **kwargs):
        from ase.calculators.siesta import Siesta
        command = '{} < PREFIX.fdf > PREFIX.out'.format(self.executable)
        return Siesta(command=command, pseudo_path=str(self.pseudo_path),
                      **kwargs)

    @classmethod
    def fromconfig(cls, config):
        paths = config.datafiles['siesta']
        assert len(paths) == 1
        path = paths[0]
        return cls(config.executables['siesta'], str(path))


class ExternalCodes:
    def __init__(self, executables, datafiles):
        self.executables = executables
        self.datafiles = datafiles
        self._factories = {}

    def __getitem__(self, name):
        # Hmm.  We could also just return a new factory instead of
        # caching them.
        if name not in self._factories:
            cls = factory_classes[name]
            try:
                factory = cls.fromconfig(self)
            except KeyError:
                pytest.skip('Missing configuration for {}'.format(name))
            self._factories[name] = factory
        return self._factories[name]


@pytest.fixture(scope='session')
def factories(executables, datafiles, enabled_calculators):
    return ExternalCodes(executables, datafiles)


def newfactory(name):
    @pytest.fixture(scope='session')
    def factory(factories):
        return factories[name]
    factory.__name__ = '{}_factory'.format(name)
    return factory


abinit_factory = newfactory('abinit')
cp2k_factory = newfactory('cp2k')
espresso_factory = newfactory('espresso')
siesta_factory = newfactory('siesta')


@pytest.fixture
def code(request, factories):
    return factories[request.param]


def pytest_generate_tests(metafunc):
    """Parametrize tests using our custom markers.

    We want tests marked with @pytest.mark.calculator(names) to be
    parametrized over the named calculator or calculators."""
    # Is it okay to use metafunc.definition.iter_markers()?
    # Or is that private or something?
    for marker in metafunc.definition.iter_markers(name='calculator'):
        calculator_names = marker.args[0]
        if isinstance(calculator_names, str):
            calculator_names = [calculator_names]
        metafunc.parametrize('code', calculator_names, indirect=True)
