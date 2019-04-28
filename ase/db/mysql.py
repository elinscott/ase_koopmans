from pymysql import connect

from ase.db import SQLite3Database
from ase.db import init_statements
from ase.db import index_statements


class Connection(object):
    def __init__(self, host=None, user=None, passwd=None,
                 db_name=None, sshtunnel=None):
        self.sshtunnel = sshtunnel
        self.con = connect(host=self.host, user=self.user, passwd=self.passwd,
                           db=self.db_name)

    def cursor(self):
        return self.con.cursor()

    def commit(self):
        self.con.commit()

    def close(self):
        self.con.close()


class MySQLDatabase(SQLite3Database):
    type = 'mysql'
    default = 'DEFAULT'

    def __init__(self, filename=None, create_indices=True,
                 use_lock_file=False, serial=False):
                super(self).__init__(filename, create_indices, use_lock_file,
                                     serial)

                self.host = None
                self.username = None
                self.passwd = None
                self.db_name = None
                self._parse_filename(filename)

    def _parse_filename(self, filename):
        filename = filename.replace('mysql://', '')

        splitted = filename.split(':')
        self.host = splitted[0]
        self.username = splitted[1]
        self.passwd = splitted[2]
        self.db_name = splitted[3]

    def _connect(self):
        return Connection(host=self.host, user=self.user, passwd=self.passwd,
                          db_name=self.db_name)

    def _initialize(self, con):
        if self.initialized:
            return

        cur = con.cursor()

        information_exists = True
        try:
            cur.execute("SELECT 1 FROM 'information'")
        except Exception as exc:
            print(type(exc))
            information_exists = False

        if not information_exists:
            # We need to initialize the DB
            sql = schema_update(init_statements)
            cur.execute(sql)

            if self.create_indices:
                cur.execute(index_statements)
            con.commit()

        self.initialized = True


def schema_update(sql):
    for a, b in [('REAL', 'DOUBLE'),
                 ('INTEGER PRIMARY KEY AUTOINCREMENT',
                  'INT NOT NULL AUTO_INCREMENT')]:
        sql = sql.replace(a, b)

    txt2jsonb = ['calculator_parameters', 'key_value_pairs', 'data']

    for column in txt2jsonb:
        sql = sql.replace('{} TEXT,'.format(column),
                          '{} JSON,'.format(column))

    return sql

