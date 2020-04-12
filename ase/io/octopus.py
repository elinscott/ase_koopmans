import os
from ase.utils import reader
from ase.calculators.octopus import parse_input_file, kwargs2atoms


@reader
def read_octopus(fd, get_kwargs=False):
    kwargs = parse_input_file(fd)

    # input files may contain internal references to other files such
    # as xyz or xsf.  We need to know the directory where the file
    # resides in order to locate those.  If fd is a real file
    # object, it contains the path and we can use it.  Else assume
    # pwd.
    #
    # Maybe this is ugly; maybe it can lead to strange bugs if someone
    # wants a non-standard file-like type.  But it's probably better than
    # failing 'ase gui somedir/inp'
    try:
        fname = fd.name
    except AttributeError:
        directory = None
    else:
        directory = os.path.split(fname)[0]

    atoms, remaining_kwargs = kwargs2atoms(kwargs, directory=directory)
    if get_kwargs:
        return atoms, remaining_kwargs
    else:
        return atoms
