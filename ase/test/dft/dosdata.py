import pytest

from ase.dft.dosdata import DOSData

# class TestDosData:
#     """Test the abstract base class for DOS data"""

sample_info = [(None, {}),
               ({}, {}),
               ({'symbol': 'C', 'index': '2', 'strangekey': 'isallowed'},
                {'symbol': 'C', 'index': '2', 'strangekey': 'isallowed'}),
               ('notadict', TypeError),
               (False, TypeError)]

@pytest.mark.parametrize('info, expected', sample_info)
def test_dd_init_info(info, expected):
    """Check 'info' parameter is handled properly"""
    if isinstance(expected, type) and isinstance(expected(), Exception):
        with pytest.raises(expected):
            dos_data = DOSData(info=info)
    else:
        dos_data = DOSData(info=info)
        assert dos_data.info == expected

