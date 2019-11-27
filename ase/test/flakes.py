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
          'E241,W291,E302,E225,E128,E261,E251,E265,E226,E231,E129,'
          'E741,W293,W503,W504,E122,W391,E126,E501')

ignore = 'W504,W291,W503,E741'

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
    'calculators/demonnano.py',
    'test/qbox/qboxdata.py',
    'test/eam_pot.py',
    'test/fio/vasp_out.py']


def flakes(path):
    # We use pyflakes because it's fast.  flake8 takes more than a minute.
    #
    # However pyflakes is not very configurable.  We want to ignore some
    # statements, namely those commented '# noqa'.
    # Hence we do some parsing.
    print('flake8:', path)
    proc = Popen([sys.executable, '-m', 'flake8', str(path),
                  # '--ignore', ignore,
                  '--exclude', ','.join(str(path / x) for x in exclude),
                  '-j', '1'],
                 stdout=PIPE)
    stdout, stderr = proc.communicate()
    stdout = stdout.decode('utf8')

    trouble_lines = []
    from collections import defaultdict
    E = defaultdict(int)
    F = defaultdict(int)
    X = {}
    for stdout_line in stdout.splitlines():
        tokens = stdout_line.split(':', 3)
        filename, lineno, colno, complaint = tokens
        lineno = int(lineno)
        e = complaint.strip().split()[0]
        E[e] += 1
        X[e] = complaint
        # if e == 'E501':
        F[filename] += 1
        """
        # Find offending line and check if it has the noqa comment:
        with open(filename) as fd:
            for i in range(lineno):
                line = next(fd)

            if re.search(r'#\s*noqa', line.lower()):
                print(filename, 'ignore', line.strip())
            else:
                print(filename, 'trouble', line.strip())
                trouble_lines.append(stdout_line)
        """
    for e, n in sorted(E.items(), key=lambda i: i[1]):
        print(e, n, X[e])
    for f, n in sorted(F.items(), key=lambda i: i[1]):
        print(f, n)
    msg = 'Flakes:\n{}'.format('\n'.join(trouble_lines))
    #assert not trouble_lines, msg


flakes(asepath)
# flakes(asedocpath)
