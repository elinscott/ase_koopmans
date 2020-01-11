import os
import pytest
from ase.utils import workdir


@pytest.fixture(scope='session', autouse=True)
def disable_calculators(request):
    from ase.test.testsuite import disable_calculators
    from ase.calculators.calculator import names as calculator_names
    enabled_names = os.environ.get('ASE_TEST_CALCULATORS', '').split()
    disable_calculators([name for name in calculator_names
                         if name not in enabled_names])


@pytest.fixture(autouse=True)
def module_workdir(request, tmp_path_factory, worker_id):
    modname = request.module.__name__
    assert modname.startswith('ase.test.')
    testdirname = modname[len('ase.test.'):]

    if testdirname == 'test_modules':
        funcname = request.function.__name__
        assert funcname.startswith('test_')
        testdirname = funcname[len('test_'):]

    basetemp = tmp_path_factory.getbasetemp()
    if worker_id != 'master':
        # pytest-xdist creates a subdirectory for each worker.
        # We would rather just have our dirs named the same no matter what.
        basetemp = basetemp.parent
    tempdir = basetemp / testdirname
    with workdir(tempdir, mkdir=True):
        yield tempdir
