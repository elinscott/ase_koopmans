from pymysql import connect
from pymysql.err import ProgrammingError
from copy import deepcopy

from ase.db.sqlite import SQLite3Database
from ase.db.sqlite import init_statements
from ase.db.sqlite import VERSION


class Connection(object):
    def __init__(self, host=None, user=None, passwd=None,
                 db_name=None):
        self.con = connect(host=host, user=user, passwd=passwd, db=db_name)

    def cursor(self):
        return MySQLCursor(self.con.cursor())

    def commit(self):
        self.con.commit()

    def close(self):
        self.con.close()


class MySQLCursor(object):
    # MySQL has some reserved words that is not reserved
    # in SQLite. Hence, we need to redefine those words
    # As these words are hardcoded into the code, we cannot
    # simply rename them during initialization
    table_redefines = {
        'keys': 'attribute_keys',
    }

    field_redefines = {
        'key': 'attribute'
    }

    invalid_mysql_tables = ['number_key_values', 'keys', 'text_key_values']

    def __init__(self, cur):
        self.cur = cur

    def _is_select_statement(self, sql):
        return sql.lower().startswith('select')

    def _is_update_statement(self, sql):
        return sql.lower().startswith('update')

    def _is_insert_statement(self, sql):
        return sql.lower().startswith('insert into')

    def _is_delete_statement(self, sql):
        return sql.lower().startswith('delete')

    def _is_known_statement(self, sql):
        return self._is_select_statement(sql) or \
            self._is_update_statement(sql) or \
            self._is_insert_statement(sql) or \
            self._is_delete_statement(sql)

    def _redefine_invalid_tables(self, sql):
        for invalid in self.invalid_mysql_tables:
            if invalid in sql:
                sql = sql.replace(
                    'key=', '{}='.format(self.field_redefines['key']))
        return sql

    def execute(self, sql, params=None):
        if ' keys ' in sql:
            if not self._is_known_statement(sql):
                raise ValueError('{} is unknown'.format(sql))
            sql = sql.replace(
                ' keys ', ' {} '.format(self.table_redefines['keys']))

        sql = sql.replace('?', '%s')
        if params is None:
            params = ()
        self.cur.execute(sql, params)

    def fetchone(self):
        return self.cur.fetchone()

    def fetchall(self):
        return self.cur.fetchall()

    def executemany(self, sql, values):
        sql = self._redefine_invalid_tables(sql)
        if ' keys ' in sql:
            if not self._is_known_statement(sql):
                raise ValueError('{} is unknown'.format(sql))
            sql = sql.replace(
                ' keys ', ' {} '.format(self.table_redefines['keys']))
        sql = sql.replace('?', '%s')
        self.cur.executemany(sql, values)


class MySQLDatabase(SQLite3Database):
    type = 'mysql'
    default = 'DEFAULT'

    def __init__(self, filename=None, create_indices=True,
                 use_lock_file=False, serial=False):
                super(MySQLDatabase, self).__init__(
                    filename, create_indices, use_lock_file, serial)

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
        return Connection(host=self.host, user=self.username,
                          passwd=self.passwd, db_name=self.db_name)

    def _initialize(self, con):
        if self.initialized:
            return

        cur = con.cursor()

        information_exists = True
        try:
            cur.execute("SELECT 1 FROM information")
        except ProgrammingError as exc:
            information_exists = False

        if not information_exists:
            # We need to initialize the DB
            # MySQL require that id is explicitly set as primary key
            # in the systems table
            init_statements[0] = init_statements[0][:-1] + ', PRIMARY KEY(id))'

            statements = schema_update(deepcopy(init_statements))
            for statement in statements:
                cur.execute(statement)

            if self.create_indices:
                print("Warning! The MySQL implementation does currently "
                      "not support indexing because of the datatype TEXT "
                      "cannot be hashed.")
            con.commit()
            self.version = VERSION
        else:
            cur.execute('select * from information')

            for name, value in cur.fetchall():
                if name == 'version':
                    self.version = int(value)

        self.initialized = True

    def blob(self, array):
        if array is None:
            return None
        return super(MySQLDatabase, self).blob(array).tobytes()

    def get_last_id(self, cur):
        cur.execute("SELECT AUTO_INCREMENT FROM information_schema.TABLES "
                    "WHERE TABLE_SCHEMA = '{}' AND TABLE_NAME = 'systems'"
                    "".format(self.db_name))
        return cur.fetchone()[0] - 1


def schema_update(statements):
    for i, statement in enumerate(statements):
        for a, b in [('REAL', 'DOUBLE'),
                     ('INTEGER PRIMARY KEY AUTOINCREMENT',
                      'INT NOT NULL AUTO_INCREMENT')]:
            statements[i] = statement.replace(a, b)

    # MySQL does not support UNIQUE constraint on TEXT
    # need to use VARCHAR. The unique_id is generated with
    # randint(16**31, 16**32-1) so it will contain 32
    # hex-characters
    statements[0] = statements[0].replace('TEXT UNIQUE', 'VARCHAR(32) UNIQUE')
    statements[2] = statements[2].replace('keys', 'attribute_keys')

    txt2jsonb = ['calculator_parameters', 'key_value_pairs', 'data']

    for column in txt2jsonb:
        statements[0] = statements[0].replace(
            '{} TEXT,'.format(column),
            '{} JSON,'.format(column))

    tab_with_key_field = ['attribute_keys', 'number_key_values',
                          'text_key_values']

    for i, statement in enumerate(statements):
        for tab in tab_with_key_field:
            if tab in statement:
                statements[i] = statement.replace(
                    'key TEXT', 'attribute_key TEXT')
    return statements
