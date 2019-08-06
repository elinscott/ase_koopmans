from io import StringIO, BytesIO
from ase.io import iread, write


def atoms_to_string(images, format=None, **kwargs):
    """Convert atoms or multiple atoms objects to string."""
    return _atoms_to_thing(images, str, format, **kwargs)


def atoms_to_bytes(images, format=None, **kwargs):
    """Convert atoms or multiple atoms objects to bytes."""
    return _atoms_to_thing(images, bytes, format, **kwargs)


def parse_images(string, format=None, **kwargs):
    """Parse string or bytes into list of atoms objects."""
    if isinstance(string, str):
        buf = StringIO(string)
    else:
        assert isinstance(string, bytes)
        buf = BytesIO(string)
    images = list(iread(buf, format=format, **kwargs))
    return images


def _atoms_to_thing(images, buftype, format, **kwargs):
    if buftype == str:
        buf = StringIO()
    elif buftype == bytes:
        assert buftype == bytes
        buf = BytesIO()
    write(buf, images, format=format, **kwargs)
    txt = buf.getvalue()
    assert isinstance(txt, buftype)
    return txt
