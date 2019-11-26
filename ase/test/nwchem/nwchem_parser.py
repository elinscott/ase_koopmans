from ase.io.nwchem.nwreader import _pattern_test_data

for regex, pattern in _pattern_test_data:
    assert regex.match(pattern) is not None, pattern
