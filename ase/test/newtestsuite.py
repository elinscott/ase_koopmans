from typing import Dict, Any, List
import runpy
import unittest
from pathlib import Path
import ase.test as asetest
from ase.utils import workdir


ignorefiles = {'__init__.py', 'testsuite.py', 'newtestsuite.py'}


def find_all_test_modules() -> List[str]:
    testdir = Path(asetest.__file__).parent
    testfiles = sorted(testdir.glob('*.py'))
    testfiles += sorted(testdir.glob('*/*.py'))

    modules = []
    for testfile in testfiles:
        if testfile.name in ignorefiles:
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


disabled_normally = {'abinit',
                     'ace', 'aims', 'aims', 'amber',
                     'calculator',
                     'calculators',
                     'castep', 'cp2k', 'crystal',
                     'demon', 'demonnano',
                     'dftb', 'dmol', 'elk', 'espresso',
                     'exciting', 'fleur',
                     'gaussian', 'gpaw', 'gromacs', 'jacapo',
                     'kim', 'lammpslib', 'lammpsrun', 'nwchem',
                     'octopus', 'onetep', 'openmx', 'psi4',
                     'qbox', 'qchem', 'siesta', 'turbomole', 'vasp'}


def define_all_tests(namespace: Dict[str, Any]):
    for module in find_all_test_modules():
        assert '-' not in module, module
        tokens = module.split('.')
        assert tokens[0] == 'ase'
        assert tokens[1] == 'test'
        if tokens[2] in disabled_normally:
            continue
        testfunc = define_script_test_function(module)
        namespace[testfunc.__name__] = testfunc


define_all_tests(globals())
