import pytest
from ase.spacegroup import (get_bravais_class,
                            get_point_group,
                            polar_space_group)
import ase.lattice


functions = [get_bravais_class, get_point_group, polar_space_group]


@pytest.mark.parametrize("sg,lattice,point_group,polar",
                         [[100, ase.lattice.TET, '4mm', True],
                          [225, ase.lattice.FCC, '4/m -3 2/m', False]])
def test_valid_spacegroup(sg, lattice, point_group, polar):
    assert get_bravais_class(sg) == lattice
    assert get_point_group(sg) == point_group
    assert polar_space_group(sg) == polar


@pytest.mark.parametrize("func", functions)
def test_nonpositive_spacegroup(func):
    with pytest.raises(ValueError, match="positive"):
        func(0)


@pytest.mark.parametrize("func", functions)
def test_bad_spacegroup(func):
    with pytest.raises(ValueError, match="Bad"):
        func(400)
