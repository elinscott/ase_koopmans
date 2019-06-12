"""Test to ensure that md logger and trajectory contain same data"""
from pathlib import Path

try:
    import pytest
except ModuleNotFoundError:
    pass

import numpy as np

from ase.optimize import FIRE, BFGS
from ase.data import s22
from ase.calculators.tip3p import TIP3P
from ase.constraints import FixBondLengths
from ase.md.verlet import VelocityVerlet
from ase.md.langevin import Langevin
import ase.units as u
from ase.io import Trajectory

md_cls_and_kwargs = [
    (VelocityVerlet, {}),
    (Langevin, {"temperature": 300 * u.kB, "friction": 0.02}),
]


# prepare atoms object for testing
dimer = s22.create_s22_system("Water_dimer")

calc = TIP3P(rc=9.0)
dimer.constraints = FixBondLengths(
    [(3 * i + j, 3 * i + (j + 1) % 3) for i in range(2) for j in [0, 1, 2]]
)


def fmax(forces):
    return np.sqrt((forces ** 2).sum(axis=1).max())


@pytest.mark.parametrize("cls", [FIRE, BFGS])
@pytest.mark.parametrize("atoms", [dimer])
@pytest.mark.parametrize("calc", [calc])
def test_opt(cls, atoms, calc, logfile="opt.log", trajectory="opt.traj"):
    """run optimization and verify that log and trajectory coincide"""

    clear_files()

    opt_atoms = atoms.copy()
    opt_atoms.constraints = atoms.constraints

    opt_atoms.calc = calc

    opt = cls(opt_atoms, logfile=logfile, trajectory=trajectory)

    # Run optimizer two times
    opt.run(0.2)
    opt.run(0.1)

    with Trajectory(trajectory) as traj, open(logfile) as f:
        next(f)
        for _, (a, line) in enumerate(zip(traj, f)):
            fmax1 = float(line.split()[-1])
            fmax2 = fmax(a.get_forces())

            assert np.allclose(fmax1, fmax2, atol=0.01), (fmax1, fmax2)


@pytest.mark.parametrize("cls, kwargs", md_cls_and_kwargs)
@pytest.mark.parametrize("atoms", [dimer])
@pytest.mark.parametrize("calc", [calc])
def test_md(
    cls, atoms, calc, kwargs, logfile="md.log", timestep=1 * u.fs, trajectory="md.traj"
):
    """ run MD for 10 steps and verify that trajectory and log coincide """

    clear_files()

    if hasattr(atoms, "constraints"):
        del atoms.constraints

    atoms.calc = calc

    md = cls(atoms, logfile=logfile, timestep=timestep, trajectory=trajectory, **kwargs)

    # run md two times
    md.run(steps=5)
    md.run(steps=5)

    # assert log file has correct length
    length = sum(1 for l in open(logfile))
    assert length == 12, length

    with Trajectory(trajectory) as traj, open(logfile) as f:
        next(f)
        for _, (a, line) in enumerate(zip(traj, f)):
            Epot1, T1 = float(line.split()[-3]), float(line.split()[-1])
            Epot2, T2 = a.get_potential_energy(), a.get_temperature()

            assert np.allclose(T1, T2, atol=0.1), (T1, T2)
            assert np.allclose(Epot1, Epot2, atol=0.01), (Epot1, Epot2)


# clear previous runs
def clear_files(suffixes=None, verbose=False):
    """clear all files with given suffix"""
    if suffixes is None:
        suffixes = ["log", "traj"]

    for suffix in suffixes:
        files = Path().glob(f"*.{suffix}")
        for file in files:
            if file.exists():
                file.unlink()
                if verbose:
                    print(f".. {file} removed.")


if __name__ == "__main__":
    clear_files()

    # test optimizer
    for cls in (FIRE, BFGS):
        test_opt(cls, dimer, calc)

    # test md
    del dimer.constraints
    for cls, kwargs in md_cls_and_kwargs:
        test_md(cls, dimer, calc, kwargs=kwargs)
