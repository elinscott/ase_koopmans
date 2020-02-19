import numpy as np
import pytest

from ase.dft.dosdata import DOSData, RawDOSData

class TestDosData:
    """Test the abstract base class for DOS data"""

    sample_info = [(None, {}),
                   ({}, {}),
                   ({'symbol': 'C', 'index': '2', 'strangekey': 'isallowed'},
                    {'symbol': 'C', 'index': '2', 'strangekey': 'isallowed'}),
                   ('notadict', TypeError),
                   (False, TypeError)]

    @pytest.mark.parametrize('info, expected', sample_info)
    def test_dosdata_init_info(self, info, expected):
        """Check 'info' parameter is handled properly"""
        if isinstance(expected, type) and isinstance(expected(), Exception):
            with pytest.raises(expected):
                dos_data = DOSData(info=info)
        else:
            dos_data = DOSData(info=info)
            assert dos_data.info == expected

    dosdata_abc_notimplemented_methods_args = [('get_energies', tuple()),
                                               ('get_weights', tuple()),
                                               ('sample', ([0.1, 0.2],)),
                                               ]
    @pytest.mark.parametrize('method, args',
                             dosdata_abc_notimplemented_methods_args)
    def test_dosdata_notimplemented(self, method, args):
        """Check NotImplementedError raised from abstract base class"""
        dos_data = DOSData()
        with pytest.raises(NotImplementedError):
            getattr(dos_data, method)(*args)

class TestRawDosData:
    """Test the raw DOS data container"""

    @pytest.fixture
    def sparse_dos(self):
        return RawDOSData([1.2, 3.4, 5.], [3., 2.1, 0.],
                          info={'symbol': 'H', 'number': '1', 'food': 'egg'})

    @pytest.fixture
    def another_sparse_dos(self):
        return RawDOSData([8., 2., 2., 5.], [1., 1., 1., 1.],
                          info={'symbol': 'H', 'number': '2'})

    def test_raw_data_access(self, sparse_dos):
        assert sparse_dos.info == {'symbol': 'H', 'number': '1', 'food': 'egg'}
        assert np.allclose(sparse_dos.get_energies(), [1.2, 3.4, 5.])
        assert np.allclose(sparse_dos.get_weights(), [3., 2.1, 0.])
        
    def test_raw_data_addition(self, sparse_dos, another_sparse_dos):
        summed_dos = sparse_dos + another_sparse_dos
        assert summed_dos.info == {'symbol': 'H'}
        assert np.allclose(summed_dos.get_energies(),
                           [1.2, 3.4, 5., 8., 2., 2., 5.])
        assert np.allclose(summed_dos.get_weights(),
                           [3., 2.1, 0., 1., 1., 1., 1.])
