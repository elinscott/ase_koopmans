import os

import pytest

from ase.db import connect
from ase.calculators.emt import EMT
from ase.build import molecule


@pytest.fixture(scope='module')
def url():
    pytest.importorskip('pymysql')
    ON_CI_SERVER = 'CI_PROJECT_DIR' in os.environ

    if ON_CI_SERVER:
        URL = 'mysql://root:ase@mysql:3306/testase_mysql'
        # HOST = 'mysql'
        # USER = 'root'
        # PASSWD = 'ase'
        # DB_NAME = 'testase_mysql'
    else:
        URL = os.environ.get('MYSQL_DB_URL')
        # HOST = os.environ.get('MYSQL_HOST', None)
        # USER = os.environ.get('MYSQL_USER', None)
        # PASSWD = os.environ.get('MYSQL_PASSWD', None)
        # DB_NAME = os.environ.get('MYSQL_DB_NAME', None)

    if URL is None:
        pytest.skip('Not on GitLab CI server. To run this test '
                    'host, username, password and database name '
                    'must be in the environment variables '
                    'MYSQL_HOST, MYSQL_USER, MYSQL_PASSWD and '
                    'MYSQL_DB_NAME, respectively.')


@pytest.fixture
def db(url):
    return connect(url)


@pytest.fixture
def h2o():
    return molecule('H2O')


def test_connect(db):
    db.delete([row.id for row in db.select()])


def test_write_read(db):
    co = molecule('CO')
    uid = db.write(co, tag=1, type='molecule')

    co_db = db.get(id=uid)
    atoms_db = co_db.toatoms()

    assert len(atoms_db) == 2
    assert atoms_db[0].symbol == co[0].symbol
    assert atoms_db[1].symbol == co[1].symbol
    assert co_db.tag == 1
    assert co_db.type == 'molecule'


def test_write_read_with_calculator(db, h2o):
    calc = EMT(dummy_param=2.4)
    h2o.set_calculator(calc)

    uid = db.write(h2o)

    h2o_db = db.get(id=uid).toatoms(attach_calculator=True)

    calc_db = h2o_db.get_calculator()
    assert calc_db.parameters['dummy_param'] == 2.4

    # Check that get_atoms function works
    db.get_atoms(H=2)


def test_update(db, h2o):
    uid = db.write(h2o, type='molecule')
    db.update(id=uid, type='oxide')

    atoms_type = db.get(id=uid).type

    assert atoms_type == 'oxide'


def test_delete(db, h2o):
    uid = db.write(h2o, type='molecule')

    # Make sure that we can get the value
    db.get(id=uid)
    db.delete([uid])

    with pytest.raises(KeyError):
        db.get(id=uid)


def test_read_write_bool_key_value_pair(db, h2o):
    # Make sure we can read and write boolean key value pairs
    uid = db.write(h2o, is_water=True, is_solid=False)
    row = db.get(id=uid)
    assert row.is_water
    assert not row.is_solid
