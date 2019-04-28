from ase.db import connect

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

test_connect()
