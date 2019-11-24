# Run pyflakes on main source dir and documentation.
import sys
from pathlib import Path
import re
from subprocess import Popen, PIPE
import unittest
import ase

try:
    import pyflakes  # noqa
except ImportError:
    raise unittest.SkipTest('pyflakes module not available')


asepath = Path(ase.__path__[0])
asedocpath = asepath.parent / 'doc'

ignore = ('E722,E303,E221,E731,E203,E262,E202,E402,E502,E127,E201,E305,'
          'E241,W291,E302,E225,E128,E261,E251,E501,E265,E226,E231,E129,'
          'E741,W293,W503,W504,E122,W391,E126')

exclude = [
    'calculators/jacapo/*',
    'data/*',
    'transport/tools.py',
    'calculators/gulp.py',
    'dimer.py',
    'io/pov.py',
    'test/combine_mm.py',
    'optimize/fmin_bfgs.py',
    'utils/linesearch.py',
    'transport/calculators.py',
    'utils/ff.py',
    'calculators/ase_qmmm_manyqm.py',
    'utils/memory.py',
    'optimize/oldqn.py',
    'calculators/demonnano.py']


def flakes(path):
    # We use pyflakes because it's fast.  flake8 takes more than a minute.
    #
    # However pyflakes is not very configurable.  We want to ignore some
    # statements, namely those commented '# noqa'.
    # Hence we do some parsing.
    print('flake8:', path)
    print([sys.executable, '-m', 'flake8', str(path),
           '--ignore', ignore,
           '--exclude', ','.join(str(path / x) for x in exclude)])
    proc = Popen([sys.executable, '-m', 'flake8', str(path),
                  '--ignore', ignore,
                  '--exclude', ','.join(str(path / x) for x in exclude)],
                 stdout=PIPE)
    stdout, stderr = proc.communicate()
    stdout = stdout.decode('utf8')

    trouble_lines = []
    for stdout_line in stdout.splitlines():
        tokens = stdout_line.split(':', 2)
        if len(tokens) != 3:
            1 / 0  # asfjhcontinue
        filename, lineno, complaint = tokens
        lineno = int(lineno)

        # Find offending line and check if it has the noqa comment:
        with open(filename) as fd:
            for i in range(lineno):
                line = next(fd)

            if re.search(r'#\s*noqa', line.lower()):
                print(filename, 'ignore', line.strip())
            else:
                print(filename, 'trouble', line.strip())
                trouble_lines.append(stdout_line)

    msg = 'Flakes:\n{}'.format('\n'.join(trouble_lines))
    assert not trouble_lines, msg


flakes(asepath)
# flakes(asedocpath)
