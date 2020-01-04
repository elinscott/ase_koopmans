from ase.test.newtestsuite import TestModule

# Define unittests for all the old-style tests which are visible
# to pytest as individual test functions:
TestModule.add_oldstyle_tests_to_namespace(globals())
