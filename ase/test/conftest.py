import os
import pytest
from ase.utils import workdir
from pathlib import Path


@pytest.fixture(scope='session', autouse=True)
def disable_calculators(request):
    from ase.test.testsuite import disable_calculators
    from ase.calculators.calculator import names as calculator_names
    enabled_names = os.environ.get('ASE_TEST_CALCULATORS', '').split()
    disable_calculators([name for name in calculator_names
                         if name not in enabled_names])


@pytest.fixture(autouse=True)
def use_tmp_workdir(tmp_path):
    # Pytest can on some systems provide a Path from pathlib2.  Normalize:
    path = Path(str(tmp_path))
    with workdir(path, mkdir=True):
        yield tmp_path
