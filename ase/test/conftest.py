import os
import pytest
from ase.utils import workdir
from pathlib import Path


@pytest.fixture(scope='session', autouse=True)
def disable_calculators(request):
    from ase.test.testsuite import disable_calculators, test_calculator_names
    from ase.calculators.calculator import names as calculator_names
    enabled_names = os.environ.get('ASE_TEST_CALCULATORS', '').split()
    test_calculator_names += enabled_names
    disable_calculators([name for name in calculator_names
                         if name not in enabled_names])


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
    asetest = pytest.importorskip('asetest')
    return asetest.datafiles


@pytest.fixture(scope='session')
def executables():
    import json
    path = os.environ.get('ASE_EXECUTABLE_CONFIGFILE')
    if path is None:
        path = Path.home() / '.ase' / 'executables.json'
    path = Path(path)
    if not path.exists():
        pytest.skip('Environment variable ASE_EXECUTABLE_CONFIG not set '
                    'and {} does not exist'.format(path))
    dct = json.loads(path.read_text())
    return dct


class AbinitFactory:
    def __init__(self, executable, pp_paths):
        self.executable = executable
        self.pp_paths = pp_paths

    def _base_kw(self):
        exe = self.executable
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


class CP2KFactory:
    def __init__(self, executable):
        self.executable = executable

    def calc(self, **kwargs):
        from ase.calculators.cp2k import CP2K
        return CP2K(command=self.executable, **kwargs)


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

@pytest.fixture(scope='session')
def abinit_factory(executables, datafiles):
    return AbinitFactory(executables['abinit'],
                         datafiles.paths['abinit'])


@pytest.fixture(scope='session')
def cp2k_factory(executables):
    return CP2KFactory(executables['cp2k'])


@pytest.fixture(scope='session')
def espresso_factory(executables, datafiles):
    paths = datafiles.paths['espresso']
    assert len(paths) == 1
    return EspressoFactory(executables['espresso'], paths[0])
