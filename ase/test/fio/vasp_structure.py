import os
import numpy as np
import numpy.testing
import ase
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

    def assert_atoms_almost_equal(self, atoms, other):
        """Compare two Atoms objects, raising AssertionError if different

        This is quite redundant with the Atoms.__eq__ method, but more
        tolerant of minor differences in the positions by a) comparing scaled
        (i.e. wrapped) positions and b) using a comparison that tolerates small
        errors.

        """
        self.assertIsInstance(atoms, ase.Atoms)
        self.assertIsInstance(other, ase.Atoms)
        self.assertEqual(len(atoms), len(other))
        self.assertEqual(atoms.arrays['numbers'].tolist(),
                         other.arrays['numbers'].tolist())
        numpy.testing.assert_array_almost_equal(
            atoms.get_scaled_positions(wrap=True),
            other.get_scaled_positions(wrap=True),
            decimal=12)
        numpy.testing.assert_array_almost_equal(atoms.cell, other.cell,
                                                decimal=12)
        self.assertTrue((atoms.pbc == other.pbc).all())

    def assert_trajectory_almost_equal(self, traj1, traj2):
        self.assertEqual(len(traj1), len(traj2))
        for image, other in zip(traj1, traj2):
            self.assert_atoms_almost_equal(image, other)

    def test_roundtrip(self):
        # Create a series of translated cells
        trajectory = [self.NaCl.copy() for i in range(5)]
        for i, atoms in enumerate(trajectory):
            atoms.set_scaled_positions(atoms.get_scaled_positions()
                                       + i * np.array([0.05, 0, 0.02]))
            atoms.wrap()

        ase.io.write(self.outfile, trajectory, format='vasp-xdatcar')
        roundtrip_trajectory = ase.io.read(self.outfile, index=':')
        self.assert_trajectory_almost_equal(trajectory, roundtrip_trajectory)

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


if __name__ in ['__main__', 'test']:
    runner = unittest.TextTestRunner(resultclass=XdatcarTestResults)
    runner.run(suite())
