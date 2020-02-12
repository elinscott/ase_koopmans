from ase.build import fcc111
from ase.ga.slab_operators import (CutSpliceSlabCrossover,
                                   RandomCompositionMutation)


def get_Cu_slab():
    a = 1
    size = (2, 4, 3)
    p1 = fcc111('Cu', size, orthogonal=True, a=a)
    p1.info['confid'] = 1
    return p1


def test_cut_splice():
    ratio = .4
    op = CutSpliceSlabCrossover(min_ratio=ratio)
    p1 = get_Cu_slab()
    natoms = len(p1)

    p2 = get_Cu_slab()
    p2.symbols = ['Au'] * natoms

    p2.info['confid'] = 2
    child, desc = op.get_new_individual([p1, p2])
    assert desc == 'CutSpliceSlabCrossover: Parents 1 2'

    # Check the ratio of elements
    syms = child.get_chemical_symbols()
    new_ratio = syms.count('Au') / natoms
    assert new_ratio > ratio and new_ratio < 1 - ratio

    op = CutSpliceSlabCrossover(element_pools=['Cu', 'Au'],
                                allowed_compositions=[(12, 12)])
    child = op.operate(p1, p2)
    assert child.get_chemical_symbols().count('Au') == 12


def test_random_mutation():
    p1 = get_Cu_slab()
    p1.symbols[3] = 'Au'
    op = RandomCompositionMutation(element_pools=['Cu', 'Au'],
                                   allowed_compositions=[(12, 12),
                                                         (18, 6)])
    child, _ = op.get_new_individual([p1])
    assert child.get_chemical_symbols().count('Au') in [6, 12]


if __name__ == "__main__":
    test_random_mutation()
