# Run flake8 tests on main source dir and documentation.
import sys
import pathlib
import unittest
from subprocess import Popen, PIPE
import ase

try:
    import flake8  # noqa
except ImportError:
    raise unittest.SkipTest('flake8 module not available')


asepath = pathlib.Path(ase.__path__[0])
asedocpath = asepath.parent / 'doc'

# API of flake8 is not stable according to the project's documentation.
# We therefore call it through a subprocess.

def flakes(path):
    print('flake8:', path)
    proc = Popen([sys.executable, '-m',
                  'flake8', '--ignore', 'E,W', str(path)],
                 stdout=PIPE)
    stdout, stderr = proc.communicate()
    code = proc.wait()
    assert code == 0, 'Flakes:\n' + stdout.decode('utf8', 'replace')


flakes(asepath)
flakes(asedocpath)
