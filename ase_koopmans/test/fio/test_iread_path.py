from pathlib import Path
import pytest
from ase_koopmans.io import iread


def test_iread_path():
    path = Path('hello')
    iterator = iread(path)
    with pytest.raises(FileNotFoundError):
        next(iterator)
