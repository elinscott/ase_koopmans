import pytest

from ase_koopmans import Atoms
from ase_koopmans.build import bulk, fcc111
from ase_koopmans.io import write


@pytest.fixture(
    params=[
        bulk('Ti'),
        fcc111('Au', size=(1, 1, 1)),
        pytest.param(Atoms('H', cell=[0, 0, 1], pbc=[0, 0, 1]),
                     marks=pytest.mark.xfail),  # Please_koopmans investigate failure
    ],
    ids=lambda atoms: f'{atoms.cell.rank}-dim',
)
def file(request):
    atoms = request.param
    file = f'atoms.{atoms.cell.rank}dim.traj'
    write(file, atoms)
    return file


def test_bzplot(cli, file, plt):
    cli.ase_koopmans(['reciprocal', file, 'bandpath.svg'])
