import numpy as np
from ase.geometry.cell import complete_cell

eps = 1E-10

def random_unit_vector():
    while 1:
        v = np.random.uniform(-1, 1, 3)
        norm = np.linalg.norm(v)
        if norm > eps:
            return v / norm

np.random.seed(0)

for it in range(100):

    cell = np.zeros((3, 3))
    index = np.random.randint(0, 3)
    cell[index] = random_unit_vector()
    complete = complete_cell(cell)

    assert abs(np.linalg.norm(complete[index]) - 1) < eps
    assert np.linalg.det(complete) > 0
    assert np.linalg.matrix_rank(complete) == 3

