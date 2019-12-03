import os
import numpy as np
import ase.build
import ase.io

import unittest

class Test_xdatcar_roundtrip(unittest.TestCase):
    def setUp(self):
        self.outfile = 'NaCl.XDATCAR'

        self.NaCl = ase.build.bulk('NaCl', 'rocksalt', a=5.64)

    def tearDown(self):
        if os.path.isfile(self.outfile):
            os.remove(self.outfile)

    def test_roundtrip(self):
        # Create a series of translated cells
        trajectory = [self.NaCl.copy() for i in range(5)]
        for i, atoms in enumerate(trajectory):
            atoms.set_scaled_positions(atoms.get_scaled_positions()
                                           + i * np.array([0.05, 0, 0.02]))
            atoms.wrap()

        ase.io.write(self.outfile, trajectory, format='vasp-xdatcar')

        roundtrip_trajectory = ase.io.read(self.outfile, index=':')

        self.assertEqual(roundtrip_trajectory, trajectory)


def suite():
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(
        Test_xdatcar_roundtrip)
    return suite

# Instead of keeping/displaying unittest results, escalate errors so ASE unit
# test system can handle them
class XdatcarTestResults(unittest.TestResult):
    def addFailure(self, test, err):
        raise err[1]

    def addError(self, test, err):
        raise err[1]

runner = unittest.TextTestRunner(resultclass=XdatcarTestResults)
runner.run(suite())
