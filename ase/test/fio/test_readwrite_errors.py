def test_readwrite_errors():
    from io import StringIO
    from ase.io import read, write
    from ase.build import bulk
    from ase.test import must_raise
    from ase.io.formats import UnknownFileTypeError

    atoms = bulk('Au')
    fd = StringIO()

    with must_raise(UnknownFileTypeError):
        write(fd, atoms, format='hello')

    with must_raise(UnknownFileTypeError):
        read(fd, format='hello')
