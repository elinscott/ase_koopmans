import pytest
from typing import Iterable

from ase.dft.doscollection import DOSCollection
from ase.dft.dosdata import DOSData, RawDOSData


class MinimalDOSCollection(DOSCollection):
    """Inherit from abstract base class to check its features"""
    def __init__(self, dos_series: Iterable[DOSData]) -> None:
        super().__init__(dos_series)


class TestDOSCollection:
    @pytest.fixture
    def rawdos(self):
        return RawDOSData([1., 2., 4.], [2., 3., 2.],
                          info={'my_key': 'my_value'})
    @pytest.fixture
    def another_rawdos(self):
        return RawDOSData([3., 2., 4.], [1., 0., 2.],
                          info={'other_key': 'other_value'})

    @pytest.mark.parametrize('n_entries', [0, 1, 3])
    def test_sequence(self, rawdos, n_entries):
        dc = DOSCollection([rawdos] * n_entries)
        assert len(dc) == n_entries
        for i in range(n_entries):
            assert dc[i] == rawdos

        with pytest.raises(IndexError):
            dc[n_entries + 1]
        with pytest.raises(TypeError):
            dc['hello']

    equality_data = [([], [], True),
                     ([rawdos], [rawdos], True),
                     ([rawdos, another_rawdos],
                      [rawdos, another_rawdos], True),
                     ([], [rawdos], False),
                     ([rawdos], [], False),
                     ([rawdos, another_rawdos], [rawdos], False),
                     ([rawdos, another_rawdos],
                      [another_rawdos, rawdos], False)]

    @pytest.mark.parametrize('series_1, series_2, isequal', equality_data)
    def test_collection_equality(self, rawdos, another_rawdos,
                                 series_1, series_2, isequal):
        assert (DOSCollection(series_1) == DOSCollection(series_2)) == isequal

    @pytest.mark.parametrize('other', [True, 1, 0.5, 'string', rawdos])
    def test_equality_wrongtype(self, rawdos, other):
        assert not DOSCollection([rawdos]) == other

    def test_addition(self, rawdos, another_rawdos):
        dc = DOSCollection([rawdos])

        double_dc = dc + dc
        assert len(double_dc) == 2
        assert double_dc[0] == rawdos
        assert double_dc[1] == rawdos

        assert dc + DOSCollection([another_rawdos]) == dc + another_rawdos
