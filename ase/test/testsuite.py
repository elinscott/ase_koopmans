import os
import sys
import subprocess
from contextlib import contextmanager
import importlib
import unittest
import warnings
import argparse

import pytest

from ase.calculators.calculator import names as calc_names, get_calculator_class
from ase.cli.info import print_info


def importorskip(module):
    # The pytest importorskip() function raises a wierd exception
    # which claims to come from the builtins module, but doesn't!
    #
    # That exception messes with our pipeline when sending stacktraces
    # through multiprocessing.  Argh.
    #
    # We provide our own implementation then!
    try:
        return importlib.import_module(module)
    except ImportError:  # From py3.6 we can use ModuleNotFoundError
        raise unittest.SkipTest('Optional module not present: {}'
                                .format(module))


test_calculator_names = ['emt']


def require(calcname):
    if calcname not in test_calculator_names:
        raise unittest.SkipTest('use --calculators={0} to enable'
                                .format(calcname))


def disable_calculators(names):
    for name in names:
        if name in ['emt', 'lj', 'eam', 'morse', 'tip3p']:
            continue
        try:
            cls = get_calculator_class(name)
        except ImportError:
            pass
        else:
            def get_mock_init(name):
                def mock_init(obj, *args, **kwargs):
                    raise unittest.SkipTest('use --calculators={0} to enable'
                                            .format(name))
                return mock_init

            def mock_del(obj):
                pass
            cls.__init__ = get_mock_init(name)
            cls.__del__ = mock_del


def cli(command, calculator_name=None):
    if (calculator_name is not None and
        calculator_name not in test_calculator_names):
        return
    actual_command = ' '.join(command.split('\n')).strip()
    proc = subprocess.Popen(actual_command,
                            shell=True,
                            stdout=subprocess.PIPE)
    print(proc.stdout.read().decode())
    proc.wait()

    if proc.returncode != 0:
        raise RuntimeError('Command "{}" exited with error code {}'
                           .format(actual_command, proc.returncode))


class must_raise:
    """Context manager for checking raising of exceptions."""
    def __init__(self, exception):
        self.exception = exception

    def __enter__(self):
        pass

    def __exit__(self, exc_type, exc_value, tb):
        if exc_type is None:
            raise RuntimeError('Failed to fail: ' + str(self.exception))
        return issubclass(exc_type, self.exception)


@contextmanager
def must_warn(category):
    with warnings.catch_warnings(record=True) as ws:
        yield
        did_warn = any(w.category == category for w in ws)
    if not did_warn:
        raise RuntimeError('Failed to warn: ' + str(category))


@contextmanager
def no_warn():
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore')
        yield


MULTIPROCESSING_MAX_WORKERS = 32
MULTIPROCESSING_DISABLED = 0
MULTIPROCESSING_AUTO = -1


class CLICommand:
    """Run ASE's test-suite.

    By default, tests for external calculators are skipped.  Enable with
    "-c name".
    """

    @staticmethod
    def add_arguments(parser):
        parser.add_argument(
            '-c', '--calculators',
            help='comma-separated list of calculators to test')
        parser.add_argument('--list', action='store_true',
                            help='print all tests and exit')
        parser.add_argument('--list-calculators', action='store_true',
                            help='print all calculator names and exit')
        parser.add_argument('-j', '--jobs', type=int, metavar='N',
                            default=MULTIPROCESSING_AUTO,
                            help='number of worker processes.  '
                            'By default use all available processors '
                            'up to a maximum of {}.  '
                            '0 disables multiprocessing'
                            .format(MULTIPROCESSING_MAX_WORKERS))
        parser.add_argument('-v', '--verbose', action='store_true',
                            help='write test outputs to stdout.  '
                            'Mostly useful when inspecting a single test')
        parser.add_argument('--strict', action='store_true',
                            help='convert warnings to errors')
        parser.add_argument('--nogui', action='store_true',
                            help='do not run graphical tests')
        parser.add_argument('tests', nargs='*',
                            help='specify particular test files.  '
                            'Glob patterns are accepted.')
        parser.add_argument('--pytest', nargs=argparse.REMAINDER,
                            help='forward all remaining arguments to pytest.  '
                            'See pytest --help')

    @staticmethod
    def run(args):
        if args.calculators:
            calculators = args.calculators.split(',')
            os.environ['ASE_TEST_CALCULATORS'] = ' '.join(calculators)
        else:
            calculators = []

        print_info()

        if args.list_calculators:
            for name in calc_names:
                print(name)
            sys.exit(0)

        for calculator in calculators:
            if calculator not in calc_names:
                sys.stderr.write('No calculator named "{}".\n'
                                 'Possible CALCULATORS are: '
                                 '{}.\n'.format(calculator,
                                                ', '.join(calc_names)))
                sys.exit(1)

        if args.nogui:
            os.environ.pop('DISPLAY')

        pytest_args = ['-v']

        def add_args(*args):
            pytest_args.extend(args)

        if args.list:
            add_args('--collect-only')

        if args.jobs == MULTIPROCESSING_DISABLED:
            pass
        elif args.jobs == MULTIPROCESSING_AUTO:
            add_args('--numprocesses=auto',
                     '--maxprocesses={}'.format(MULTIPROCESSING_MAX_WORKERS))
        else:
            add_args('--numprocesses={}'.format(args.jobs))

        add_args('--pyargs')

        if args.tests:
            from ase.test.newtestsuite import TestModule

            dct = TestModule.all_test_modules_as_dict()

            for testname in args.tests:
                mod = dct[testname]
                if mod.is_pytest_style:
                    pytest_args.append(mod.module)
                else:
                    # XXX Not totally logical
                    add_args('ase.test.test_modules::{}'
                             .format(mod.pytest_function_name))
        else:
            add_args('ase.test')

        if args.verbose:
            add_args('--capture=no')

        if args.pytest:
            add_args(*args.pytest)

        print()
        calcstring = ','.join(calculators) if calculators else 'none'
        print('Enabled calculators: {}'.format(calcstring))
        print()
        print('About to run pytest with these parameters:')
        for line in pytest_args:
            print('    ' + line)
        exitcode = pytest.main(pytest_args)
        sys.exit(exitcode)
