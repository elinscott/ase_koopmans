"""Old test suite expressed as pytest-type functions.

This module exposes the old test suite as module-level test_xxx()
functions recognized by pytest.

We should start writing unittests as functions (or classes) that
pytest works with.  We can also port the old tests to that form, where
it matters.

Either way: The goal is that the list of tests provided by this module
becomes shorter, and the list of tests in other modules becomes
longer."""


from typing import Dict, Any, Generator
import runpy
from pathlib import Path
import ase.test as asetest
from ase.utils import workdir


from ase.calculators.calculator import names as calculator_names

# Ignore calculator tests (for now):
#calculators = {
#    'abinit', 'ace', 'aims', 'aims', 'amber',
#    'calculator', 'calculators',
#    'castep', 'cp2k', 'crystal', 'demon', 'demonnano',
#    'dftb', 'dmol', 'elk', 'espresso',
#    'exciting', 'fleur', 'gaussian', 'gpaw', 'gromacs', 'jacapo',
#    'kim', 'lammpslib', 'lammpsrun', 'nwchem',
#    'octopus', 'onetep', 'openmx', 'psi4',
#    'qbox', 'qchem', 'siesta', 'turbomole', 'vasp',
#}


class TestModule:
    ignorefiles = {'__init__.py', 'testsuite.py', 'newtestsuite.py',
                   'conftest.py'}
    testdir = Path(asetest.__file__).parent

    def __init__(self, testname: str):
        # Testname is e.g. "fio.dftb".
        # Referring to the file ase/test/fio/dftb.py
        # or the module ase.test.fio.dftb
        self.testname = testname

    @property
    def module(self) -> str:
        return 'ase.test.{}'.format(self.testname)

    @property
    def path(self) -> Path:
        return Path(self.module.replace('.', '/') + '.py')

    @property
    def is_pytest_style(self) -> bool:
        # Files named test_* or *_test are are picked up by pytest
        # automatically.  We call these "pytest-style" modules.
        #
        # The other modules must be for our own old test suite.
        name = self.path.name
        return 'test_' in name or '_test' in name

    @classmethod
    def from_relpath(cls, relpath) -> "TestModule":
        module = TestModule.relpath2module(relpath)
        return TestModule(module)

    @staticmethod
    def filename_to_testname(filename):
        name, py = str(filename).rsplit('.', 1)
        assert py == 'py'
        return name.replace('/', '.')

    def __repr__(self):
        return 'TestModule({})'.format(self.module)

    @classmethod
    def glob_all_test_modules(cls) -> "Generator[TestModule]":
        """Return a list of modules ['ase.test.xxx', 'ase.test.yyy', ...]."""
        testfiles = sorted(cls.testdir.glob('*.py'))
        testfiles += sorted(cls.testdir.glob('*/*.py'))
        # XXX Some tests were added at */*/*.py level, but the old test suite
        # never globbed so deep.  So these tests never ran.
        # We can/should rehabilitate them.

        for testfile in testfiles:
            if testfile.name in cls.ignorefiles:
                continue
            #if testfile.parent.name in calculator_names:
            #    continue
            #if testfile.parent.name in ['calculators', 'calculator']:
            #    continue
            if '#' in testfile.name:
                continue  # Ignore certain backup files.
            rel_testfile = testfile.relative_to(cls.testdir)
            testname = cls.filename_to_testname(rel_testfile)
            yield TestModule(testname)

    def define_script_test_function(self):
        module = self.module
        pytestname = 'test_' + self.testname.replace('.', '_')

        def test_script(tmp_path):
            with workdir(tmp_path):
                runpy.run_module(module, run_name='test')

        test_script.__name__ = pytestname
        return test_script

    @classmethod
    def add_oldstyle_tests_to_namespace(cls, namespace: Dict[str, Any]):
        for testmodule in cls.glob_all_test_modules():
            if testmodule.is_pytest_style:
                continue

            testfunc = testmodule.define_script_test_function()
            namespace[testfunc.__name__] = testfunc
