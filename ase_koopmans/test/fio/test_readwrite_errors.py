def test_readwrite_errors():
    from io import StringIO
    from ase_koopmans.io import read, write
    from ase_koopmans.build import bulk
    from ase_koopmans.test import must_raise
    from ase_koopmans.io.formats import UnknownFileTypeError

    atoms = bulk('Au')
    fd = StringIO()

    with must_raise(UnknownFileTypeError):
        write(fd, atoms, format='hello')

    with must_raise(UnknownFileTypeError):
        read(fd, format='hello')
