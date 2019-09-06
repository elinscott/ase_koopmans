"""File formats.

This module implements the read(), iread() and write() functions in ase.io.
For each file format there is a namedtuple (IOFormat) that has the following
elements:

* a read(filename, index, **kwargs) generator that will yield Atoms objects
* a write(filename, images) function
* a 'single' boolean (False if multiple configurations is supported)
* a 'acceptsfd' boolean (True if file-descriptors are accepted)

There is a dict 'ioformats' that is filled with IOFormat objects as they are
needed.  The 'initialize()' function will create the IOFormat object by
looking at the all_formats dict and by importing the correct read/write
functions from the correct module.  The 'single' and 'acceptsfd' bools are
parsed from two-charcter string in the all_formats dict below.


Example
=======

The xyz format is implemented in the ase/io/xyz.py file which has a
read_xyz() generator and a write_xyz() function.

"""

import functools
import inspect
import os
import sys

from ase.atoms import Atoms
from ase.utils import import_module, basestring, PurePath, lazyproperty
from ase.parallel import parallel_function, parallel_generator


class UnknownFileTypeError(Exception):
    pass


class IOFormat:
    def __init__(self, name, desc, code, module_name, extensions):
        self.name = name
        self.description = desc
        assert len(code) == 2
        assert code[0] in list('+1')
        assert code[1] in list('BFS')
        self.code = code
        self.module_name = module_name
        self.extensions = extensions

    def __getitem__(self, i):
        return (self.description, self.code)[i]

    @property
    def single(self):
        return self.code[0] == '1'

    @property
    def _formatname(self):
        return self.name.replace('-', '_')

    @property
    def read(self):
        read = getattr(self.module, 'read_' + self._formatname, None)
        if read and not inspect.isgeneratorfunction(read):
            read = functools.partial(wrap_read_function, read)
        return read

    @property
    def write(self):
        return getattr(self.module, 'write_' + self._formatname, None)

    @property
    def acceptsfd(self):
        return self.code[1] != 'S'

    @property
    def isbinary(self):
        return self.code[1] == 'B'

    @lazyproperty
    def module(self):
        try:
            return import_module(self.module_name)
        except ImportError as err:
            raise UnknownFileTypeError('File format not recognized: %s.  '
                                       'Error: %s' % (format, err))


ioformats = {}  # will be filled at run-time

# 1=single, +=multiple, F=accepts a file-descriptor, S=needs a file-name str,
# B=like F, but opens in binary mode
all_formats = ioformats  # XXX We should keep one of these.


extension2format = {}

def define_format(name, desc, code, *, module=None, extensions=None):
    if module is None:
        module = name.replace('-', '_')
    if extensions is None:
        extensions = []
    elif isinstance(extensions, str):
        extensions = [extensions]

    fmt = IOFormat(name, desc, code, module_name='ase.io.' + module,
                   extensions=extensions)

    for ext in extensions:
        if ext in extension2format:
            raise ValueError('extension "{}" already registered'.format(ext))
        extension2format[ext] = fmt

    ioformats[name] = fmt
    return fmt

F = define_format
F('abinit', 'ABINIT input file', '1F'),
F('aims', 'FHI-aims geometry file', '1S'),
F('aims-output', 'FHI-aims output', '+S',
  module='aims', extensions='in'),
F('bundletrajectory', 'ASE bundle trajectory', '+S'),
F('castep-castep', 'CASTEP output file', '+F',
  module='castep', extensions='castep'),
F('castep-cell', 'CASTEP geom file', '1F',
  module='castep', extensions='cell'),
F('castep-geom', 'CASTEP trajectory file', '+F',
  module='castep', extensions='geom'),
F('castep-md', 'CASTEP molecular dynamics file', '+F',
  module='castep', extensions='md'),
F('castep-phonon', 'CASTEP phonon file', '1F',
  module='castep', extensions='phonon'),
F('cfg', 'AtomEye configuration', '1F'),
F('cif', 'CIF-file', '+B'),
F('cmdft', 'CMDFT-file', '1F'),
F('cp2k-dcd', 'CP2K DCD file', '+B',
  module='cp2k', extensions='dcd'),
F('crystal', 'Crystal fort.34 format', '1S',
  extensions=['f34', '34']),
F('cube', 'CUBE file', '1F'),
F('dacapo', 'Dacapo netCDF output file', '1F'),
F('dacapo-text', 'Dacapo text output', '1F',
  module='dacapo'),
F('db', 'ASE SQLite database file', '+S'),
F('dftb', 'DftbPlus input file', '1S'),
F('dlp4', 'DL_POLY_4 CONFIG file', '1F',
  module='dlp4', extensions='config'),
F('dlp-history', 'DL_POLY HISTORY file', '+F',
  module='dlp4'),
F('dmol-arc', 'DMol3 arc file', '+S',
  module='dmol'),
F('dmol-car', 'DMol3 structure file', '1S',
  module='dmol', extensions='car'),
F('dmol-incoor', 'DMol3 structure file', '1S',
  module='dmol'),
F('elk', 'ELK atoms definition', '1S'),
F('eon', 'EON CON file', '+F',
  extensions='con'),
F('eps', 'Encapsulated Postscript', '1S'),
F('espresso-in', 'Quantum espresso in file', '1F',
  module='espresso', extensions='pwi'),
F('espresso-out', 'Quantum espresso out file', '+F',
  module='espresso', extensions=['out', 'pwo']),
F('etsf', 'ETSF format', '1S'),
F('exciting', 'exciting input', '1S',
  extensions='exi'),
F('extxyz', 'Extended XYZ file', '+F'),
F('findsym', 'FINDSYM-format', '+F'),
F('gaussian', 'Gaussian com (input) file', '1S',
  extensions='com'),
F('gaussian-out', 'Gaussian output file', '1F',
  module='gaussian', extensions='log'),
F('acemolecule-out', 'ACE output file', '1S',
  module='acemolecule'),
F('acemolecule-input', 'ACE input file', '1S',
  module='acemolecule'),
F('gen', 'DFTBPlus GEN format', '1F'),
F('gif', 'Graphics interchange format', '+S',
  module='animation'),
F('gpaw-out', 'GPAW text output', '+F'),
F('gpw', 'GPAW restart-file', '1S'),
F('gromacs', 'Gromacs coordinates', '1S',
  extensions='gro'),
F('gromos', 'Gromos96 geometry file', '1F', extensions='g96'),
F('html', 'X3DOM HTML', '1F', module='x3d'),
F('iwm', '?', '1F'),
F('json', 'ASE JSON database file', '+F', module='db'),
F('jsv', 'JSV file format', '1F'),
F('lammps-dump', 'LAMMPS dump file', '+F', module='lammpsrun'),
F('lammps-data', 'LAMMPS data file', '1F', module='lammpsdata'),
F('magres', 'MAGRES ab initio NMR data file', '1F'),
F('mol', 'MDL Molfile', '1F'),
F('mp4', 'MP4 animation', '+S',
  module='animation'),
F('mustem', 'muSTEM xtl file', '1F',
  extensions='stl'),
F('mysql', 'ASE MySQL database file', '+S',
  module='db'),
F('netcdftrajectory', 'AMBER NetCDF trajectory file', '+S'),
F('nomad-json', 'JSON from Nomad archive', '+F'),
F('nwchem', 'NWChem input file', '1F',
  extensions='nw'),
F('octopus', 'Octopus input file', '1F'),
F('proteindatabank', 'Protein Data Bank', '+F',
  extensions='pdb'),
F('png', 'Portable Network Graphics', '1S'),
F('postgresql', 'ASE PostgreSQL database file', '+S', module='db'),
F('pov', 'Persistance of Vision', '1S'),
F('py', 'Python file', '+F'),
F('qbox', 'QBOX output file', '+F'),
F('res', 'SHELX format', '1S', extensions='shelx'),
F('rmc6f', 'RMCProfile', '1S'),
F('sdf', 'SDF format', '1F'),
F('struct', 'WIEN2k structure file', '1S', module='wien2k'),
F('struct_out', 'SIESTA STRUCT file', '1F', module='siesta'),
F('traj', 'ASE trajectory', '+B', module='trajectory'),
F('trj', 'Old ASE pickle trajectory', '+S',
  module='pickletrajectory'),
F('turbomole', 'TURBOMOLE coord file', '1F'),
F('turbomole-gradient', 'TURBOMOLE gradient file', '+F',
  module='turbomole'),
F('v-sim', 'V_Sim ascii file', '1F', extensions='ascii'),
F('vasp', 'VASP POSCAR/CONTCAR file', '1F',
  extensions='poscar'),
F('vasp-out', 'VASP OUTCAR file', '+F', module='vasp'),
F('vasp-xdatcar', 'VASP XDATCAR file', '+F', module='vasp'),
F('vasp-xml', 'VASP vasprun.xml file', '+F', module='vasp'),
F('vti', 'VTK XML Image Data', '1F', module='vtkxml'),
F('vtu', 'VTK XML Unstructured Grid', '1F', module='vtkxml'),
F('x3d', 'X3D', '1S'),
F('xsd', 'Materials Studio file', '1F'),
F('xsf', 'XCrySDen Structure File', '+F'),
F('xtd', 'Materials Studio file', '+F'),
F('xyz', 'XYZ-file', '+F')


netcdfconventions2format = {
    'http://www.etsf.eu/fileformats': 'etsf',
    'AMBER': 'netcdftrajectory'
}


def get_ioformat(format):
    """Initialize and return IOFormat tuple."""
    #initialize(format)
    return ioformats[format]


def get_compression(filename):
    """
    Parse any expected file compression from the extension of a filename.
    Return the filename without the extension, and the extension. Recognises
    ``.gz``, ``.bz2``, ``.xz``.

    >>> get_compression('H2O.pdb.gz')
    ('H2O.pdb', 'gz')
    >>> get_compression('crystal.cif')
    ('crystal.cif', None)

    Parameters
    ==========
    filename: str
        Full filename including extension.

    Returns
    =======
    (root, extension): (str, str or None)
        Filename split into root without extension, and the extension
        indicating compression format. Will not split if compression
        is not recognised.
    """
    # Update if anything is added
    valid_compression = ['gz', 'bz2', 'xz']

    # Use stdlib as it handles most edge cases
    root, compression = os.path.splitext(filename)

    # extension keeps the '.' so remember to remove it
    if compression.strip('.') in valid_compression:
        return root, compression.strip('.')
    else:
        return filename, None


def open_with_compression(filename, mode='r'):
    """
    Wrapper around builtin `open` that will guess compression of a file
    from the filename and open it for reading or writing as if it were
    a standard file.

    Implemented for ``gz``(gzip), ``bz2``(bzip2) and ``xz``(lzma).

    Supported modes are:
       * 'r', 'rt', 'w', 'wt' for text mode read and write.
       * 'rb, 'wb' for binary read and write.

    Parameters
    ==========
    filename: str
        Path to the file to open, including any extensions that indicate
        the compression used.
    mode: str
        Mode to open the file, same as for builtin ``open``, e.g 'r', 'w'.

    Returns
    =======
    fd: file
        File-like object open with the specified mode.
    """

    # Compressed formats sometimes default to binary, so force
    # text mode:
    if mode == 'r':
        mode = 'rt'
    elif mode == 'w':
        mode = 'wt'
    elif mode == 'a':
        mode = 'at'

    root, compression = get_compression(filename)

    if compression is None:
        return open(filename, mode)
    elif compression == 'gz':
        import gzip
        fd = gzip.open(filename, mode=mode)
    elif compression == 'bz2':
        import bz2
        fd = bz2.open(filename, mode=mode)
    elif compression == 'xz':
        try:
            from lzma import open as lzma_open
        except ImportError:
            from backports.lzma import open as lzma_open
        fd = lzma_open(filename, mode)
    else:
        fd = open(filename, mode)

    return fd


def wrap_read_function(read, filename, index=None, **kwargs):
    """Convert read-function to generator."""
    if index is None:
        yield read(filename, **kwargs)
    else:
        for atoms in read(filename, index, **kwargs):
            yield atoms


def write(filename, images, format=None, parallel=True, append=False,
          **kwargs):
    """Write Atoms object(s) to file.

    filename: str or file
        Name of the file to write to or a file descriptor.  The name '-'
        means standard output.
    images: Atoms object or list of Atoms objects
        A single Atoms object or a list of Atoms objects.
    format: str
        Used to specify the file-format.  If not given, the
        file-format will be taken from suffix of the filename.
    parallel: bool
        Default is to write on master only.  Use parallel=False to write
        from all slaves.
    append: bool
        Default is to open files in 'w' or 'wb' mode, overwriting
        existing files.  In some cases opening the file in 'a' or 'ab'
        mode (appending) is usefull,
        e.g. writing trajectories or saving multiple Atoms objects in one file.
        WARNING: If the file format does not support multiple entries without
        additional keywords/headers, files created using 'append=True'
        might not be readable by any program! They will nevertheless be
        written without error message.

    The use of additional keywords is format specific."""

    if isinstance(filename, PurePath):
        filename = str(filename)

    if isinstance(filename, basestring):
        filename = os.path.expanduser(filename)
        fd = None
        if filename == '-':
            fd = sys.stdout
            filename = None
        elif format is None:
            format = filetype(filename, read=False)
    else:
        fd = filename
        filename = None

    format = format or 'json'  # default is json

    io = get_ioformat(format)

    _write(filename, fd, format, io, images, parallel=parallel, append=append,
           **kwargs)


@parallel_function
def _write(filename, fd, format, io, images, parallel=None, append=False,
           **kwargs):
    if isinstance(images, Atoms):
        images = [images]

    if io.single:
        if len(images) > 1:
            raise ValueError('{}-format can only store 1 Atoms object.'
                             .format(format))
        images = images[0]

    if io.write is None:
        raise ValueError("Can't write to {}-format".format(format))

    # Special case for json-format:
    if format == 'json' and (len(images) > 1 or append):
        if filename is not None:
            io.write(filename, images, append=append, **kwargs)
            return
        raise ValueError("Can't write more than one image to file-descriptor "
                         'using json-format.')

    if io.acceptsfd:
        open_new = (fd is None)
        if open_new:
            mode = 'wb' if io.isbinary else 'w'
            if append:
                mode = mode.replace('w', 'a')
            fd = open_with_compression(filename, mode)
        io.write(fd, images, **kwargs)
        if open_new:
            fd.close()
    else:
        if fd is not None:
            raise ValueError("Can't write {}-format to file-descriptor"
                             .format(format))
        if 'append' in io.write.__code__.co_varnames:
            io.write(filename, images, append=append, **kwargs)
        elif append:
            raise ValueError("Cannot append to {}-format, write-function "
                             "does not support the append keyword."
                             .format(format))
        else:
            io.write(filename, images, **kwargs)


def read(filename, index=None, format=None, parallel=True, **kwargs):
    """Read Atoms object(s) from file.

    filename: str or file
        Name of the file to read from or a file descriptor.
    index: int, slice or str
        The last configuration will be returned by default.  Examples:

            * ``index=0``: first configuration
            * ``index=-2``: second to last
            * ``index=':'`` or ``index=slice(None)``: all
            * ``index='-3:`` or ``index=slice(-3, None)``: three last
            * ``index='::2`` or ``index=slice(0, None, 2)``: even
            * ``index='1::2`` or ``index=slice(1, None, 2)``: odd
    format: str
        Used to specify the file-format.  If not given, the
        file-format will be guessed by the *filetype* function.
    parallel: bool
        Default is to read on master and broadcast to slaves.  Use
        parallel=False to read on all slaves.

    Many formats allow on open file-like object to be passed instead
    of ``filename``. In this case the format cannot be auto-decected,
    so the ``format`` argument should be explicitly given."""

    if isinstance(filename, PurePath):
        filename = str(filename)
    if filename == '-':
        filename = sys.stdin
    if isinstance(index, basestring):
        try:
            index = string2index(index)
        except ValueError:
            pass

    filename, index = parse_filename(filename, index)
    if index is None:
        index = -1
    format = format or filetype(filename)
    io = get_ioformat(format)
    if isinstance(index, (slice, basestring)):
        return list(_iread(filename, index, format, io, parallel=parallel,
                           **kwargs))
    else:
        return next(_iread(filename, slice(index, None), format, io,
                           parallel=parallel, **kwargs))


def iread(filename, index=None, format=None, parallel=True, **kwargs):
    """Iterator for reading Atoms objects from file.

    Works as the `read` function, but yields one Atoms object at a time
    instead of all at once."""

    if isinstance(index, basestring):
        index = string2index(index)

    filename, index = parse_filename(filename, index)

    if index is None or index == ':':
        index = slice(None, None, None)

    if not isinstance(index, (slice, basestring)):
        index = slice(index, (index + 1) or None)

    format = format or filetype(filename)
    io = get_ioformat(format)

    for atoms in _iread(filename, index, format, io, parallel=parallel,
                        **kwargs):
        yield atoms


@parallel_generator
def _iread(filename, index, format, io, parallel=None, full_output=False,
           **kwargs):
    if isinstance(filename, basestring):
        filename = os.path.expanduser(filename)

    if not io.read:
        raise ValueError("Can't read from {}-format".format(format))

    if io.single:
        start = index.start
        assert start is None or start == 0 or start == -1
        args = ()
    else:
        args = (index,)

    must_close_fd = False
    if isinstance(filename, basestring):
        if io.acceptsfd:
            mode = 'rb' if io.isbinary else 'r'
            fd = open_with_compression(filename, mode)
            must_close_fd = True
        else:
            fd = filename
    else:
        assert io.acceptsfd
        fd = filename

    # Make sure fd is closed in case loop doesn't finish:
    try:
        for dct in io.read(fd, *args, **kwargs):
            if not isinstance(dct, dict):
                dct = {'atoms': dct}
            if full_output:
                yield dct
            else:
                yield dct['atoms']
    finally:
        if must_close_fd:
            fd.close()


def parse_filename(filename, index=None):
    if not isinstance(filename, basestring):
        return filename, index

    extension = os.path.basename(filename)
    if '@' not in extension:
        return filename, index

    newindex = None
    newfilename, newindex = filename.rsplit('@', 1)

    if isinstance(index, slice):
        return newfilename, index
    try:
        newindex = string2index(newindex)
    except ValueError:
        pass

    return newfilename, newindex


def string2index(string):
    if ':' not in string:
        return int(string)
    i = []
    for s in string.split(':'):
        if s == '':
            i.append(None)
        else:
            i.append(int(s))
    i += (3 - len(i)) * [None]
    return slice(*i)


def filetype(filename, read=True, guess=True):
    """Try to guess the type of the file.

    First, special signatures in the filename will be checked for.  If that
    does not identify the file type, then the first 2000 bytes of the file
    will be read and analysed.  Turn off this second part by using
    read=False.

    Can be used from the command-line also::

        $ ase info filename ...
    """

    ext = None
    if isinstance(filename, basestring):
        if os.path.isdir(filename):
            if os.path.basename(os.path.normpath(filename)) == 'states':
                return 'eon'
            return 'bundletrajectory'

        if filename.startswith('postgres'):
            return 'postgresql'

        if filename.startswith('mysql') or filename.startswith('mariadb'):
            return 'mysql'

        # strip any compression extensions that can be read
        root, compression = get_compression(filename)
        basename = os.path.basename(root)

        if basename == 'inp':
            return 'octopus'

        if basename.endswith('.nomad.json'):
            return 'nomad-json'

        if '.' in basename:
            ext = os.path.splitext(basename)[1].strip('.').lower()
            if ext in ['xyz', 'cube', 'json', 'cif']:
                return ext

        if 'POSCAR' in basename or 'CONTCAR' in basename:
            return 'vasp'
        if 'OUTCAR' in basename:
            return 'vasp-out'
        if 'XDATCAR' in basename:
            return 'vasp-xdatcar'
        if 'vasp' in basename and basename.endswith('.xml'):
            return 'vasp-xml'
        if basename == 'coord':
            return 'turbomole'
        if basename == 'f34':
            return 'crystal'
        if basename == '34':
            return 'crystal'
        if basename == 'gradient':
            return 'turbomole-gradient'
        if basename.endswith('I_info'):
            return 'cmdft'
        if basename == 'atoms.dat':
            return 'iwm'
        if 'CONFIG' in basename:
            return 'dlp4'
        if basename == 'HISTORY':
            return 'dlp-history'

        if not read:
            if ext is None:
                raise UnknownFileTypeError('Could not guess file type')
            return extension2format.get(ext, ext)

        fd = open_with_compression(filename, 'rb')
    else:
        fd = filename
        if fd is sys.stdin:
            return 'json'

    data = fd.read(50000)
    if fd is not filename:
        fd.close()
    else:
        fd.seek(0)

    if len(data) == 0:
        raise UnknownFileTypeError('Empty file: ' + filename)

    if data.startswith(b'CDF'):
        # We can only recognize these if we actually have the netCDF4 module.
        try:
            import netCDF4
        except ImportError:
            pass
        else:
            nc = netCDF4.Dataset(filename)
            if 'Conventions' in nc.ncattrs():
                if nc.Conventions in netcdfconventions2format:
                    return netcdfconventions2format[nc.Conventions]
                else:
                    raise UnknownFileTypeError(
                        "Unsupported NetCDF convention: "
                        "'{}'".format(nc.Conventions))
            else:
                raise UnknownFileTypeError("NetCDF file does not have a "
                                           "'Conventions' attribute.")

    for format, magic in [('traj', b'- of UlmASE-Trajectory'),
                          ('traj', b'AFFormatASE-Trajectory'),
                          ('gpw', b'- of UlmGPAW'),
                          ('gpw', b'AFFormatGPAW'),
                          ('trj', b'PickleTrajectory'),
                          ('turbomole', b'$coord'),
                          ('turbomole-gradient', b'$grad'),
                          ('dftb', b'Geometry')]:
        if data.startswith(magic):
            return format

    for format, magic in [('gpaw-out', b'  ___ ___ ___ _ _ _'),
                          ('espresso-in', b'\n&system'),
                          ('espresso-in', b'\n&SYSTEM'),
                          ('espresso-out', b'Program PWSCF'),
                          ('aims-output', b'Invoking FHI-aims ...'),
                          ('lammps-dump', b'\nITEM: TIMESTEP\n'),
                          ('qbox', b':simulation xmlns:'),
                          ('xsf', b'\nANIMSTEPS'),
                          ('xsf', b'\nCRYSTAL'),
                          ('xsf', b'\nSLAB'),
                          ('xsf', b'\nPOLYMER'),
                          ('xsf', b'\nMOLECULE'),
                          ('xsf', b'\nATOMS'),
                          ('dacapo-text',
                           b'&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n')]:
        if magic in data:
            return format

    format = extension2format.get(ext)
    if format is None and guess:
        format = ext
    if format is None:
        # Do quick xyz check:
        lines = data.splitlines()
        if lines and lines[0].strip().isdigit():
            return 'xyz'

        raise UnknownFileTypeError('Could not guess file type')

    return format


def index2range(index, nsteps):
    """Method to convert a user given *index* option to a list of indices.

    Returns a range.
    """
    if isinstance(index, int):
        if index < 0:
            tmpsnp = nsteps + index
            trbl = range(tmpsnp, tmpsnp + 1, 1)
        else:
            trbl = range(index, index + 1, 1)
    elif isinstance(index, slice):
        start = index.start
        stop = index.stop
        step = index.step

        if start is None:
            start = 0
        elif start < 0:
            start = nsteps + start

        if step is None:
            step = 1

        if stop is None:
            stop = nsteps
        elif stop < 0:
            stop = nsteps + stop

        trbl = range(start, stop, step)
    else:
        raise RuntimeError("index2range handles integers and slices only.")
    return trbl
