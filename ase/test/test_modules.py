# We rename the class so pytest does not confuse it with an actual test:
from ase.test.newtestsuite import TestModule as TstModule

# This module will appear to pytest as a long list of test functions:
TstModule.add_oldstyle_tests_to_namespace(globals())
