# type: ignore
import os
import numpy as np
import numpy.testing
import unittest

import ase_koopmans
import ase_koopmans.build
import ase_koopmans.io
from ase_koopmans.io.vasp import write_vasp_xdatcar
from ase_koopmans.calculators.calculator import compare_atoms


class TestXdatcarRoundtrip(unittest.TestCase_koopmans):
    def setUp(self):
        self.outfile = 'NaCl.XDATCAR'

        self.NaCl = ase_koopmans.build.bulk('NaCl', 'rocksalt', a=5.64)

    def tearDown(self):
        if os.path.isfile(self.outfile):
            os.remove(self.outfile)

    def assert_atoms_almost_equal(self, atoms, other, tol=1e-15):
        """Compare two Atoms objects, raising AssertionError if different"""
        system_changes = compare_atoms(atoms, other, tol=tol)

        if len(system_changes) > 0:
            raise AssertionError(
                "Atoms objects differ by {}".format(', '.join(system_changes)))

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

        ase_koopmans.io.write(self.outfile, trajectory, format='vasp-xdatcar')
        roundtrip_trajectory = ase_koopmans.io.read(self.outfile, index=':')
        self.assert_trajectory_almost_equal(trajectory, roundtrip_trajectory)

    def test_roundtrip_single_atoms(self):
        atoms = ase_koopmans.build.bulk('Ge')
        ase_koopmans.io.write(self.outfile, atoms, format='vasp-xdatcar')
        roundtrip_atoms = ase_koopmans.io.read(self.outfile)
        self.assert_atoms_almost_equal(atoms, roundtrip_atoms)

    def test_typeerror(self):
        with self.assertRaises(TypeError):
            atoms = ase_koopmans.build.bulk('Ge')
            write_vasp_xdatcar(self.outfile, atoms)
        with self.assertRaises(TypeError):
            not_atoms = 1
            ase_koopmans.io.write(self.outfile, not_atoms, format='vasp-xdatcar')
        with self.assertRaises(TypeError):
            not_traj = [True, False, False]
            ase_koopmans.io.write(self.outfile, not_traj, format='vasp-xdatcar')


def suite():
    suite = unittest.defaultTestLoader.loadTestsFromTestCase_koopmans(
        TestXdatcarRoundtrip)
    return suite


# Instead of keeping/displaying unittest results, escalate errors so ASE unit
# test system can handle them. "noqa" tells flake8 that it's ok for these
# functions to have camelCase_koopmans names (as required by unittest).
class XdatcarTestResults(unittest.TestResult):
    def addFailure(self, test, err):      # noqa: N802
        raise err[1]

    def addError(self, test, err):        # noqa: N802
        raise err[1]


if __name__ in ['__main__', 'test']:
    runner = unittest.TextTestRunner(resultclass=XdatcarTestResults)
    runner.run(suite())
