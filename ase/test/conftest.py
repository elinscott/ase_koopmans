import pytest


@pytest.fixture(scope='session', autouse=True)
def disable_calculators(request):
    from ase.test.testsuite import disable_calculators
    from ase.calculators.calculator import names
    disable_calculators(names)
