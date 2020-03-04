import pytest
import numpy as np
from ase.io.zmatrix import parse_zmatrix

tests = dict(
    h2=dict(
        zmat=('h;'
              'h 1 0.74'),
        pos=[[0, 0, 0], [0.74, 0, 0]],
    ),
    h2o2=dict(
        zmat=[
            'H1',
            'O2 1 0.9',
            'O3 2 1.4 1 105',
            'H4 3 0.9 2 105 1 120',
        ],
        pos=[[0, 0, 0],
             [0.9, 0, 0],
             [1.262, 1.352, 0],
             [1.742, 1.465, 0.753]],
    ),
    h2o_gauss=dict(
        zmat='\n'.join([
            'O',
            'H 1 rOH',
            'H 1 rOH 2 aH2O',
        ]),
        defs="rOH 0.97; aH2O 104",
        pos=[[0, 0, 0],
             [0.97, 0, 0],
             [-0.235, 0.941, 0.]],
    ),
)

tests['h2o_qchem'] = dict(
    zmat='\n'.join([
        'O1',
        'H1 O1 oh',
        'H2 O1 oh H1 hoh',
    ]),
    defs='\n'.join([
        'oh = 0.97',
        'hoh = 104.0',
    ]),
    pos=tests['h2o_gauss']['pos'],
)


@pytest.mark.parametrize('name', tests)
def test_zmatrix(name):
    atoms = parse_zmatrix(tests[name]['zmat'],
                          defs=tests[name].get('defs'))
    assert np.max(np.abs(atoms.positions - tests[name]['pos'])) < 1e-3
