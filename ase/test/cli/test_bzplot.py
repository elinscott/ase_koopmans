from subprocess import Popen

import pytest

from ase import Atoms
from ase.build import bulk, fcc111
from ase.io import write


class CLI:
    def __call__(self, args):
        proc = Popen(['ase', '-T'] + args)
        status = proc.wait()
        assert status == 0


@pytest.fixture
def cli():
    return CLI()


@pytest.fixture(
    params=[
        bulk('Ti'),
        fcc111('Au', size=(1, 1, 1)),
        pytest.param(Atoms('H', cell=[0, 0, 1], pbc=[0, 0, 1]),
                     marks=pytest.mark.xfail),  # Please investigate failure
    ],
    ids=lambda atoms: f'{atoms.cell.rank}-dim',
)
def file(request):
    atoms = request.param
    file = f'atoms.{atoms.cell.rank}dim.traj'
    write(file, atoms)
    return file


def test_bzplot(cli, file):
    cli(['reciprocal', file, 'bandpath.svg'])
