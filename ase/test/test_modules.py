"""Old test suite expressed as pytest-type functions.

This module exposes the old test suite as module-level test_xxx()
functions recognized by pytest.

We should start writing unittests as functions (or classes) that
pytest works with.  We can also port the old tests to that form, where
it matters.

Either way: The goal is that the list of tests provided by this module
becomes shorter, and the list of tests in other modules becomes
longer."""


from typing import Dict, Any, List
import runpy
from pathlib import Path
import ase.test as asetest
from ase.utils import workdir


ignorefiles = {'__init__.py', 'testsuite.py'}


# Ignore calculator tests (for now):
ignoredirs = {
    'abinit', 'ace', 'aims', 'aims', 'amber',
    'calculator', 'calculators',
    'castep', 'cp2k', 'crystal', 'demon', 'demonnano',
    'dftb', 'dmol', 'elk', 'espresso',
    'exciting', 'fleur', 'gaussian', 'gpaw', 'gromacs',
    'kim', 'lammpslib', 'lammpsrun', 'nwchem',
    'octopus', 'onetep', 'openmx', 'psi4',
    'qbox', 'qchem', 'siesta', 'turbomole', 'vasp',
}


def find_all_test_modules() -> List[str]:
    """Return a list of modules ['ase.test.xxx', 'ase.test.yyy', ...]."""
    testdir = Path(asetest.__file__).parent
    testfiles = sorted(testdir.glob('*.py'))
    testfiles += sorted(testdir.glob('*/*.py'))
    # XXX Some tests were added at */*/*.py level, but the old test suite
    # never globbed so deep.  So these tests never ran.
    # We can/should rehabilitate them.

    modules = []
    for testfile in testfiles:
        if testfile.name in ignorefiles:
            continue
        if testfile.parent.name in ignoredirs:
            continue
        if '#' in testfile.name:
            continue  # Ignore certain backup files.
        if 'test_' in testfile.name or '_test' in testfile.name:
            # These files are picked up by pytest automatically.
            # This heuristic is not perfect.
            # But we think of it as follows:
            # a test which is *not* named like this is old-style,
            # and a test named like this is pytest-style.
            continue
        rel_testfile = testfile.relative_to(testdir)
        module = str(rel_testfile).rsplit('.', 1)[0].replace('/', '.')
        modules.append('ase.test.{}'.format(module))
    return modules


def define_script_test_function(module: str):
    assert module.startswith('ase.test.')
    testname = module.split('.', 1)[1].replace('.', '_')

    def test_script(tmp_path):
        with workdir(tmp_path):
            runpy.run_module(module, run_name='test')

    test_script.__name__ = testname
    return test_script


def add_all_tests_to_namespace(namespace: Dict[str, Any]):
    for module in find_all_test_modules():
        assert '-' not in module, module  # Rename/avoid improper module names
        testfunc = define_script_test_function(module)
        namespace[testfunc.__name__] = testfunc


add_all_tests_to_namespace(globals())
