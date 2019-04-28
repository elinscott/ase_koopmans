from ase.db import connect
from ase import Atoms
from ase.calculators.emt import EMT
from ase.build import molecule

HOST = 'localhost'
USER = 'davidkl'
PASSWD = 'davidsfedora'
DB_NAME = 'demodb'


def full_db_name():
    return 'mysql://{}:{}:{}:{}'.format(HOST, USER, PASSWD, DB_NAME)


def test_connect():
    db = connect(full_db_name())

    assert db.host == HOST
    assert db.username == USER
    assert db.passwd == PASSWD
    assert db.db_name == DB_NAME


def test_write_read():
    db = connect(full_db_name())

    co = Atoms('CO', positions=[(0, 0, 0), (0, 0, 1.1)])
    uid = db.write(co, tag=1, type='molecule')

    co_db = db.get(id=uid)
    atoms_db = co_db.toatoms()

    assert len(atoms_db) == 2
    assert atoms_db[0].symbol == co[0].symbol
    assert atoms_db[1].symbol == co[1].symbol
    assert co_db.tag == 1
    assert co_db.type == 'molecule'


def test_write_read_with_calculator():
    db = connect(full_db_name())

    h2o = molecule('H2O')
    calc = EMT(dummy_param=2.4)
    h2o.set_calculator(calc)

    uid = db.write(h2o)

    h2o_db = db.get(id=uid).toatoms(attach_calculator=True)

    calc_db = h2o_db.get_calculator()
    assert calc_db.parameters['dummy_param'] == 2.4


test_connect()
test_write_read()
test_write_read_with_calculator()