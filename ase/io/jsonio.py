import datetime
import json

import numpy as np
from ase.utils import reader, writer


class MyEncoder(json.JSONEncoder):
    def default(self, obj):
        if hasattr(obj, 'todict'):
            d = obj.todict()

            if not isinstance(d, dict):
                raise RuntimeError('todict() of {} returned object of type {} '
                                   'but should have returned dict'
                                   .format(obj, type(d)))
            if hasattr(obj, 'ase_objtype'):
                d['__ase_objtype__'] = obj.ase_objtype

            return d
        if isinstance(obj, np.ndarray):
            d = {'__ase_objtype__': 'ndarray',
                 'dtype': obj.dtype.name,
                 'shape': list(obj.shape)}
            if obj.dtype == complex:
                # Cast into float array of twice the size:
                flatbuf = obj.ravel()
                flatbuf.dtype = float
                assert 2 * obj.size == flatbuf.size
                d['data'] = flatbuf.tolist()
            else:
                d['data'] = obj.ravel().tolist()
            return d
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.bool_):
            return bool(obj)
        if isinstance(obj, datetime.datetime):
            return {'__datetime__': obj.isoformat()}
        return json.JSONEncoder.default(self, obj)


encode = MyEncoder().encode


def object_hook(dct):
    if '__datetime__' in dct:
        return datetime.datetime.strptime(dct['__datetime__'],
                                          '%Y-%m-%dT%H:%M:%S.%f')
    if '__complex_ndarray__' in dct:
        r, i = (np.array(x) for x in dct['__complex_ndarray__'])
        return r + i * 1j

    if '__ase_objtype__' in dct:
        objtype = dct.pop('__ase_objtype__')

        # We just try each object type one after another and instantiate
        # them manually, depending on which kind it is.
        # We can formalize this later if it ever becomes necessary.
        if objtype == 'ndarray':
            array = np.empty(dct['shape'], dtype=dct['dtype'])
            if array.dtype == complex:
                flatbuf = array.ravel()
                flatbuf.dtype = float
                flatbuf[:] = dct['data']
            else:
                array.flat[:] = dct['data']
            return array

        elif objtype == 'cell':
            from ase.geometry.cell import Cell
            obj = Cell.new(dct['array'], dct['pbc'])
        elif objtype == 'bandstructure':
            from ase.dft.band_structure import BandStructure
            obj = BandStructure(**dct)
        elif objtype == 'bandpath':
            from ase.dft.kpoints import BandPath
            obj = BandPath(**dct)
        else:
            raise RuntimeError('Do not know how to decode object type {} '
                               'into an actual object'.format(objtype))

        assert obj.ase_objtype == objtype
        return obj

    return dct


mydecode = json.JSONDecoder(object_hook=object_hook).decode


def intkey(key):
    try:
        return int(key)
    except ValueError:
        return key


def numpyfy(obj):
    if isinstance(obj, dict):
        if '__complex_ndarray__' in obj:
            r, i = (np.array(x) for x in obj['__complex_ndarray__'])
            return r + i * 1j
        return dict((intkey(key), numpyfy(value))
                    for key, value in obj.items())
    if isinstance(obj, list) and len(obj) > 0:
        try:
            a = np.array(obj)
        except ValueError:
            pass
        else:
            if a.dtype in [bool, int, float]:
                return a
        obj = [numpyfy(value) for value in obj]
    return obj


def decode(txt, always_array=True):
    obj = mydecode(txt)
    if always_array:
        obj = numpyfy(obj)
    return obj


@reader
def read_json(fd, always_array=True):
    dct = decode(fd.read(), always_array=always_array)
    return dct


@writer
def write_json(fd, obj):
    fd.write(encode(obj))
