import functools
import json
import numbers
import operator
import os
import re
import warnings
from time import time
from typing import List, Any

import numpy as np

from ase_koopmans.atoms import Atoms
from ase_koopmans.calculators.calculator import all_properties, all_changes
from ase_koopmans.data import atomic_numbers
from ase_koopmans.db.row import AtomsRow
from ase_koopmans.formula import Formula
from ase_koopmans.io.jsonio import create_ase_koopmans_object
from ase_koopmans.parallel import world, DummyMPI, parallel_function, parallel_generator
from ase_koopmans.utils import Lock, PurePath


T2000 = 946681200.0  # January 1. 2000
YEAR = 31557600.0  # 365.25 days


# Format of key description: ('short', 'long', 'unit')
default_key_descriptions = {
    'id': ('ID', 'Uniqe row ID', ''),
    'age': ('Age', 'Time since creation', ''),
    'formula': ('Formula', 'Chemical formula', ''),
    'pbc': ('PBC', 'Periodic boundary conditions', ''),
    'user': ('Username', '', ''),
    'calculator': ('Calculator', 'ASE-calculator name', ''),
    'energy': ('Energy', 'Total energy', 'eV'),
    'natoms': ('Number of atoms', '', ''),
    'fmax': ('Maximum force', '', 'eV/Ang'),
    'smax': ('Maximum stress', 'Maximum stress on unit cell',
             '`\\text{eV/Ang}^3`'),
    'charge': ('Charge', 'Net charge in unit cell', '|e|'),
    'mass': ('Mass', 'Sum of atomic masses in unit cell', 'au'),
    'magmom': ('Magnetic moment', '', 'au'),
    'unique_id': ('Unique ID', 'Random (unique) ID', ''),
    'volume': ('Volume', 'Volume of unit cell', '`\\text{Ang}^3`')}


def now():
    """Return time since January 1. 2000 in years."""
    return (time() - T2000) / YEAR


seconds = {'s': 1,
           'm': 60,
           'h': 3600,
           'd': 86400,
           'w': 604800,
           'M': 2629800,
           'y': YEAR}

longwords = {'s': 'second',
             'm': 'minute',
             'h': 'hour',
             'd': 'day',
             'w': 'week',
             'M': 'month',
             'y': 'year'}

ops = {'<': operator.lt,
       '<=': operator.le,
       '=': operator.eq,
       '>=': operator.ge,
       '>': operator.gt,
       '!=': operator.ne}

invop = {'<': '>=', '<=': '>', '>=': '<', '>': '<=', '=': '!=', '!=': '='}

word = re.compile('[_a-zA-Z][_0-9a-zA-Z]*$')

reserved_keys = set(all_properties +
                    all_changes +
                    list(atomic_numbers) +
                    ['id', 'unique_id', 'ctime', 'mtime', 'user',
                     'fmax', 'smax',
                     'momenta', 'constraints', 'natoms', 'formula', 'age',
                     'calculator', 'calculator_parameters',
                     'key_value_pairs', 'data'])

numeric_keys = set(['id', 'energy', 'magmom', 'charge', 'natoms'])


def check(key_value_pairs):
    for key, value in key_value_pairs.items():
        if key == "external_tables":
            # Checks for external_tables are not
            # performed
            continue

        if not word.match(key) or key in reserved_keys:
            raise ValueError('Bad key: {}'.format(key))
        try:
            Formula(key, strict=True)
        except ValueError:
            pass
        else:
            warnings.warn(
                'It is best not to use keys ({0}) that are also a '
                'chemical formula.  If you do a "db.select({0!r})",'
                'you will not find rows with your key.  Instead, you wil get '
                'rows containing the atoms in the formula!'.format(key))
        if not isinstance(value, (numbers.Real, str, np.bool_)):
            raise ValueError('Bad value for {!r}: {}'.format(key, value))
        if isinstance(value, str):
            for t in [int, float]:
                if str_represents(value, t):
                    raise ValueError(
                        'Value ' + value + ' is put in as string ' +
                        'but can be interpreted as ' +
                        '{}! Please_koopmans convert '.format(t.__name__) +
                        'to {} using '.format(t.__name__) +
                        '{}(value) before '.format(t.__name__) +
                        'writing to the database_koopmans OR change ' +
                        'to a different string.')


def str_represents(value, t=int):
    try:
        t(value)
    except ValueError:
        return False
    return True


def connect(name, type='extract_from_name', create_indices=True,
            use_lock_file=True, append=True, serial=False):
    """Create connection to database_koopmans.

    name: str
        Filename or address of database_koopmans.
    type: str
        One of 'json', 'db', 'postgresql',
        (JSON, SQLite, PostgreSQL).
        Default is 'extract_from_name', which will guess the type
        from the name.
    use_lock_file: bool
        You can turn this off if you know what you are doing ...
    append: bool
        Use append=False to start a new database_koopmans.
    """

    if isinstance(name, PurePath):
        name = str(name)

    if type == 'extract_from_name':
        if name is None:
            type = None
        elif not isinstance(name, str):
            type = 'json'
        elif (name.startswith('postgresql://') or
              name.startswith('postgres://')):
            type = 'postgresql'
        elif name.startswith('mysql://') or name.startswith('mariadb://'):
            type = 'mysql'
        else:
            type = os.path.splitext(name)[1][1:]
            if type == '':
                raise ValueError('No file extension or database_koopmans type given')

    if type is None:
        return Database_koopmans()

    if not append and world.rank == 0:
        if isinstance(name, str) and os.path.isfile(name):
            os.remove(name)

    if type not in ['postgresql', 'mysql'] and isinstance(name, str):
        name = os.path.abspath(name)

    if type == 'json':
        from ase_koopmans.db.jsondb import JSONDatabase_koopmans
        return JSONDatabase_koopmans(name, use_lock_file=use_lock_file, serial=serial)
    if type == 'db':
        from ase_koopmans.db.sqlite import SQLite3Database_koopmans
        return SQLite3Database_koopmans(name, create_indices, use_lock_file,
                               serial=serial)
    if type == 'postgresql':
        from ase_koopmans.db.postgresql import PostgreSQLDatabase_koopmans
        return PostgreSQLDatabase_koopmans(name)

    if type == 'mysql':
        from ase_koopmans.db.mysql import MySQLDatabase_koopmans
        return MySQLDatabase_koopmans(name)
    raise ValueError('Unknown database_koopmans type: ' + type)


def lock(method):
    """Decorator for using a lock-file."""
    @functools.wraps(method)
    def new_method(self, *args, **kwargs):
        if self.lock is None:
            return method(self, *args, **kwargs)
        else:
            with self.lock:
                return method(self, *args, **kwargs)
    return new_method


def convert_str_to_int_float_or_str(value):
    """Safe eval()"""
    try:
        return int(value)
    except ValueError:
        try:
            value = float(value)
        except ValueError:
            value = {'True': True, 'False': False}.get(value, value)
        return value


def parse_selection(selection, **kwargs):
    if selection is None or selection == '':
        expressions = []
    elif isinstance(selection, int):
        expressions = [('id', '=', selection)]
    elif isinstance(selection, list):
        expressions = selection
    else:
        expressions = [w.strip() for w in selection.split(',')]
    keys = []
    comparisons = []
    for expression in expressions:
        if isinstance(expression, (list, tuple)):
            comparisons.append(expression)
            continue
        if expression.count('<') == 2:
            value, expression = expression.split('<', 1)
            if expression[0] == '=':
                op = '>='
                expression = expression[1:]
            else:
                op = '>'
            key = expression.split('<', 1)[0]
            comparisons.append((key, op, value))
        for op in ['!=', '<=', '>=', '<', '>', '=']:
            if op in expression:
                break
        else:
            if expression in atomic_numbers:
                comparisons.append((expression, '>', 0))
            else:
                try:
                    count = Formula(expression).count()
                except ValueError:
                    keys.append(expression)
                else:
                    comparisons.extend((symbol, '>', n - 1)
                                       for symbol, n in count.items())
            continue
        key, value = expression.split(op)
        comparisons.append((key, op, value))

    cmps = []
    for key, value in kwargs.items():
        comparisons.append((key, '=', value))

    for key, op, value in comparisons:
        if key == 'age':
            key = 'ctime'
            op = invop[op]
            value = now() - time_string_to_float(value)
        elif key == 'formula':
            if op != '=':
                raise ValueError('Use fomula=...')
            f = Formula(value)
            count = f.count()
            cmps.extend((atomic_numbers[symbol], '=', n)
                        for symbol, n in count.items())
            key = 'natoms'
            value = len(f)
        elif key in atomic_numbers:
            key = atomic_numbers[key]
            value = int(value)
        elif isinstance(value, str):
            value = convert_str_to_int_float_or_str(value)
        if key in numeric_keys and not isinstance(value, (int, float)):
            msg = 'Wrong type for "{}{}{}" - must be a number'
            raise ValueError(msg.format(key, op, value))
        cmps.append((key, op, value))

    return keys, cmps


class Database_koopmans:
    """Base_koopmans class for all database_koopmanss."""
    def __init__(self, filename=None, create_indices=True,
                 use_lock_file=False, serial=False):
        """Database_koopmans object.

        serial: bool
            Let someone else handle parallelization.  Default behavior is
            to interact with the database_koopmans on the master only and then
            distribute results to all slaves.
        """
        if isinstance(filename, str):
            filename = os.path.expanduser(filename)
        self.filename = filename
        self.create_indices = create_indices
        if use_lock_file and isinstance(filename, str):
            self.lock = Lock(filename + '.lock', world=DummyMPI())
        else:
            self.lock = None
        self.serial = serial
        self._metadata = None  # decription of columns and other stuff

    @parallel_function
    @lock
    def write(self, atoms, key_value_pairs={}, data={}, id=None, **kwargs):
        """Write atoms to database_koopmans with key-value pairs.

        atoms: Atoms object
            Write atomic numbers, positions, unit cell and boundary
            conditions.  If a calculator is attached, write also already
            calculated properties such as the energy and forces.
        key_value_pairs: dict
            Dictionary of key-value pairs.  Values must be strings or numbers.
        data: dict
            Extra stuff (not for searching).
        id: int
            Overwrite existing row.

        Key-value pairs can also be set using keyword arguments::

            connection.write(atoms, name='ABC', frequency=42.0)

        Returns integer id of the new row.
        """

        if atoms is None:
            atoms = Atoms()

        kvp = dict(key_value_pairs)  # modify a copy
        kvp.update(kwargs)

        id = self._write(atoms, kvp, data, id)
        return id

    def _write(self, atoms, key_value_pairs, data, id=None):
        check(key_value_pairs)
        return 1

    @parallel_function
    @lock
    def reserve(self, **key_value_pairs):
        """Write empty row if not already present.

        Usage::

            id = conn.reserve(key1=value1, key2=value2, ...)

        Write an empty row with the given key-value pairs and
        return the integer id.  If such a row already exists, don't write
        anything and return None.
        """

        for dct in self._select([],
                                [(key, '=', value)
                                 for key, value in key_value_pairs.items()]):
            return None

        atoms = Atoms()

        calc_name = key_value_pairs.pop('calculator', None)

        if calc_name:
            # Allow use of calculator key
            assert calc_name.lower() == calc_name

            # Fake calculator class:
            class Fake:
                name = calc_name

                def todict(self):
                    return {}

                def check_state(self, atoms):
                    return ['positions']

            atoms.calc = Fake()

        id = self._write(atoms, key_value_pairs, {}, None)

        return id

    def __delitem__(self, id):
        self.delete([id])

    def get_atoms(self, selection=None, attach_calculator=False,
                  add_additional_information=False, **kwargs):
        """Get Atoms object.

        selection: int, str or list
            See the select() method.
        attach_calculator: bool
            Attach calculator object to Atoms object (default value is
            False).
        add_additional_information: bool
            Put key-value pairs and data into Atoms.info dictionary.

        In addition, one can use keyword arguments to select specific
        key-value pairs.
        """

        row = self.get(selection, **kwargs)
        return row.toatoms(attach_calculator, add_additional_information)

    def __getitem__(self, selection):
        return self.get(selection)

    def get(self, selection=None, **kwargs):
        """Select a single row and return it as a dictionary.

        selection: int, str or list
            See the select() method.
        """
        rows = list(self.select(selection, limit=2, **kwargs))
        if not rows:
            raise KeyError('no match')
        assert len(rows) == 1, 'more than one row matched'
        return rows[0]

    @parallel_generator
    def select(self, selection=None, filter=None, explain=False,
               verbosity=1, limit=None, offset=0, sort=None,
               include_data=True, columns='all', **kwargs):
        """Select rows.

        Return AtomsRow iterator with results.  Selection is done
        using key-value pairs and the special keys:

            formula, age, user, calculator, natoms, energy, magmom
            and/or charge.

        selection: int, str or list
            Can be:

            * an integer id
            * a string like 'key=value', where '=' can also be one of
              '<=', '<', '>', '>=' or '!='.
            * a string like 'key'
            * comma separated strings like 'key1<value1,key2=value2,key'
            * list of strings or tuples: [('charge', '=', 1)].
        filter: function
            A function that takes as input a row and returns True or False.
        explain: bool
            Explain query plan.
        verbosity: int
            Possible values: 0, 1 or 2.
        limit: int or None
            Limit selection.
        offset: int
            Offset into selected rows.
        sort: str
            Sort rows after key.  Prepend with minus sign for a decending sort.
        include_data: bool
            Use include_data=False to skip reading data from rows.
        columns: 'all' or list of str
            Specify which columns from the SQL table to include.
            For example, if only the row id and the energy is needed,
            queries can be speeded up by setting columns=['id', 'energy'].
        """

        if sort:
            if sort == 'age':
                sort = '-ctime'
            elif sort == '-age':
                sort = 'ctime'
            elif sort.lstrip('-') == 'user':
                sort += 'name'

        keys, cmps = parse_selection(selection, **kwargs)
        for row in self._select(keys, cmps, explain=explain,
                                verbosity=verbosity,
                                limit=limit, offset=offset, sort=sort,
                                include_data=include_data,
                                columns=columns):
            if filter is None or filter(row):
                yield row

    def count(self, selection=None, **kwargs):
        """Count rows.

        See the select() method for the selection syntax.  Use db.count() or
        len(db) to count all rows.
        """
        n = 0
        for row in self.select(selection, **kwargs):
            n += 1
        return n

    def __len__(self):
        return self.count()

    @parallel_function
    @lock
    def update(self, id, atoms=None, delete_keys=[], data=None,
               **add_key_value_pairs):
        """Update and/or delete key-value pairs of row(s).

        id: int
            ID of row to update.
        atoms: Atoms object
            Optionally update the Atoms data (positions, cell, ...).
        data: dict
            Data dict to be added to the existing data.
        delete_keys: list of str
            Keys to remove.

        Use keyword arguments to add new key-value pairs.

        Returns number of key-value pairs added and removed.
        """

        if not isinstance(id, numbers.Integral):
            if isinstance(id, list):
                err = ('First argument must be an int and not a list.\n'
                       'Do something like this instead:\n\n'
                       'with db:\n'
                       '    for id in ids:\n'
                       '        db.update(id, ...)')
                raise ValueError(err)
            raise TypeError('id must be an int')

        check(add_key_value_pairs)

        row = self._get_row(id)
        kvp = row.key_value_pairs

        n = len(kvp)
        for key in delete_keys:
            kvp.pop(key, None)
        n -= len(kvp)
        m = -len(kvp)
        kvp.update(add_key_value_pairs)
        m += len(kvp)

        moredata = data
        data = row.get('data', {})
        if moredata:
            data.update(moredata)
        if not data:
            data = None

        if atoms:
            oldrow = row
            row = AtomsRow(atoms)
            # Copy over data, kvp, ctime, user and id
            row._data = oldrow._data
            row.__dict__.update(kvp)
            row._keys = list(kvp)
            row.ctime = oldrow.ctime
            row.user = oldrow.user
            row.id = id

        if atoms or os.path.splitext(self.filename)[1] == '.json':
            self._write(row, kvp, data, row.id)
        else:
            self._update(row.id, kvp, data)
        return m, n

    def delete(self, ids):
        """Delete rows."""
        raise NotImplementedError


def time_string_to_float(s):
    if isinstance(s, (float, int)):
        return s
    s = s.replace(' ', '')
    if '+' in s:
        return sum(time_string_to_float(x) for x in s.split('+'))
    if s[-2].isalpha() and s[-1] == 's':
        s = s[:-1]
    i = 1
    while s[i].isdigit():
        i += 1
    return seconds[s[i:]] * int(s[:i]) / YEAR


def float_to_time_string(t, long=False):
    t *= YEAR
    for s in 'yMwdhms':
        x = t / seconds[s]
        if x > 5:
            break
    if long:
        return '{:.3f} {}s'.format(x, longwords[s])
    else:
        return '{:.0f}{}'.format(round(x), s)


def object_to_bytes(obj: Any) -> bytes:
    """Serialize Python object to bytes."""
    parts = [b'12345678']
    obj = o2b(obj, parts)
    offset = sum(len(part) for part in parts)
    x = np.array(offset, np.int64)
    if not np.little_endian:
        x.byteswap(True)
    parts[0] = x.tobytes()
    parts.append(json.dumps(obj, separators=(',', ':')).encode())
    return b''.join(parts)


def bytes_to_object(b: bytes) -> Any:
    """Deserialize bytes to Python object."""
    x = np.frombuffer(b[:8], np.int64)
    if not np.little_endian:
        x = x.byteswap()
    offset = x.item()
    obj = json.loads(b[offset:].decode())
    return b2o(obj, b)


def o2b(obj: Any, parts: List[bytes]):
    if isinstance(obj, (int, float, bool, str, type(None))):
        return obj
    if isinstance(obj, dict):
        return {key: o2b(value, parts) for key, value in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [o2b(value, parts) for value in obj]
    if isinstance(obj, np.ndarray):
        assert obj.dtype != object, \
            'Cannot convert ndarray of type "object" to bytes.'
        offset = sum(len(part) for part in parts)
        if not np.little_endian:
            obj = obj.byteswap()
        parts.append(obj.tobytes())
        return {'__ndarray__': [obj.shape,
                                obj.dtype.name,
                                offset]}
    if isinstance(obj, complex):
        return {'__complex__': [obj.real, obj.imag]}
    objtype = getattr(obj, 'ase_koopmans_objtype')
    if objtype:
        dct = o2b(obj.todict(), parts)
        dct['__ase_koopmans_objtype__'] = objtype
        return dct
    raise ValueError('Objects of type {type} not allowed'
                     .format(type=type(obj)))


def b2o(obj: Any, b: bytes) -> Any:
    if isinstance(obj, (int, float, bool, str, type(None))):
        return obj

    if isinstance(obj, list):
        return [b2o(value, b) for value in obj]

    assert isinstance(obj, dict)

    x = obj.get('__complex__')
    if x is not None:
        return complex(*x)

    x = obj.get('__ndarray__')
    if x is not None:
        shape, name, offset = x
        dtype = np.dtype(name)
        size = dtype.itemsize * np.prod(shape).astype(int)
        a = np.frombuffer(b[offset:offset + size], dtype)
        a.shape = shape
        if not np.little_endian:
            a = a.byteswap()
        return a

    dct = {key: b2o(value, b) for key, value in obj.items()}
    objtype = dct.pop('__ase_koopmans_objtype__', None)
    if objtype is None:
        return dct
    return create_ase_koopmans_object(objtype, dct)
