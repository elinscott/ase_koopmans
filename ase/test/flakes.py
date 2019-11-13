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


def flakes(path):
    # We use pyflakes because it's fast.  flake8 takes more than a minute.
    #
    # However pyflakes is not very configurable.  We want to ignore some
    # statements, namely those commented '# noqa'.
    # Hence we do some parsing.
    print('pyflakes:', path)
    proc = Popen([sys.executable, '-m', 'pyflakes', str(path)],
                 stdout=PIPE)
    stdout, stderr = proc.communicate()
    stdout = stdout.decode('utf8')

    trouble_lines = []
    for stdout_line in stdout.splitlines():
        tokens = stdout_line.split(':', 2)
        if len(tokens) != 3:
            continue
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
flakes(asedocpath)
