def test_particle_comparators(seed):
    import numpy as np
    from ase_koopmans.cluster import Icosahedron
    from ase_koopmans.ga.particle_comparator import NNMatComparator
    from ase_koopmans.ga.utilities import get_nnmat
    from ase_koopmans.ga.particle_mutations import RandomPermutation

    # set up the random number generator
    rng = np.random.RandomState(seed)

    ico1 = Icosahedron('Cu', 4)
    ico1.info['confid'] = 1
    ico2 = Icosahedron('Ni', 4)
    ico1.numbers[:55] = [28] * 55
    ico2.numbers[:92] = [29] * 92

    ico1.info['data'] = {}
    ico1.info['data']['nnmat'] = get_nnmat(ico1)
    ico2.info['data'] = {}
    ico2.info['data']['nnmat'] = get_nnmat(ico2)
    comp = NNMatComparator()
    assert not comp.looks_like(ico1, ico2)

    op = RandomPermutation(rng=rng)
    a3, desc = op.get_new_individual([ico1])

    assert a3.get_chemical_formula() == ico1.get_chemical_formula()

    hard_comp = NNMatComparator(d=100)
    assert hard_comp.looks_like(ico1, a3)

    soft_comp = NNMatComparator(d=.0001)
    assert not soft_comp.looks_like(ico1, a3)
