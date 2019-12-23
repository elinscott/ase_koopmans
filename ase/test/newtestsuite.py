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


def define_all_tests(namespace: Dict[str, Any]):
    for module in find_all_test_modules():
        assert '-' not in module, module
        testfunc = define_script_test_function(module)
        namespace[testfunc.__name__] = testfunc


define_all_tests(globals())
