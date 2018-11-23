from ase.db.sqlite import SQLite3Database
from ase.parallel import parallel_function
from ase.db.core import lock, Database

class FlexibleSqlite(SQLite3Database):
    def _create_table_if_not_exists(self, name, entries):
        """Huge data structures needs to be stored differently.
           basically similar form as key_value_pairs
        """
        con = self.connection or self._connect()
        cur = con.cursor()
        dtype = self._guess_type(entries)
        sql = "CREATE TABLE IF NOT EXISTS {} (key TEXT, value {}, id INTEGER)".format(name, dtype)
        cur.execute(sql)

        # Create an index on the id for fast lookup
        sql = "CREATE INDEX IF NOT EXISTS {}_index ON {} (id)".format(name, name)
        cur.execute(sql)
        if self.connection is None:
            con.commit()
            con.close()

    def _guess_type(self, entries):
        """Guess the type based on the first entry."""
        val = entries[list(entries.keys())[0]]

        if isinstance(val, int):
            return "INTEGER"
        if isinstance(val, float):
            return "REAL"
        if isinstance(val, str):
            return "TEXT"
        raise ValueError("Unknown datatype!")

    @parallel_function
    @lock
    def insert(self, name=None, entries=None):
        """SQLite has limitations on the number of columns
           so very big tables needs a different scheme."""
        if name is None or entries is None:
            # There is nothing to do
            return
        self._create_table_if_not_exists(name, entries)
        id = entries.pop("id")

        con = self.connection or self._connect()
        cur = con.cursor()

        # First we check if entries alrady exists
        cur.execute("SELECT key FROM {} WHERE id=?".format(name), (id,))
        for item in cur:
            value = entries.pop(item[0], None)
            if value is not None:
                sql = "UPDATE {} SET value=? WHERE id=? AND key=?".format(name)
                cur.execute(sql, (value, id, item[0]))

        inserts = [(k, v, id) for k, v in entries.items()]
        sql = "INSERT INTO {} VALUES (?, ?, ?)".format(name)
        cur.executemany(sql, inserts)

        if self.connection is None:
            con.commit()
            con.close()

    def _id_exists_in_table(self, cursor, name, uid):
        """Check if the system ID already exists."""
        cursor.execute("SELECT id FROM {} " 
                       "WHERE id=? LIMIT 1".format(name), (uid,))
        return cursor.fetchone() is not None
        

    def write(self, atoms, key_value_pairs={}, data={}, id=None, **kwargs):
        extra_tables = kwargs.pop("tables", {})
        extra_tables.update(key_value_pairs.pop("tables", {}))
        uid = Database.write(self, atoms, key_value_pairs=key_value_pairs, data=data, id=id, **kwargs)
        self._insert_external_tables(extra_tables, uid)

    def _insert_external_tables(self, tables, uid):
        """Insert the external tables to the database."""
        for tabname, columns in tables.items():
            columns["id"] = uid
            self.insert(name=tabname, entries=columns)

    def update(self, id, atoms=None, delete_keys=[], data=None,
                **add_key_value_pairs):
        extra_tables = add_key_value_pairs.pop("tables", {})
        Database.update(self, id, atoms=atoms, delete_keys=delete_keys, 
                        data=data, **add_key_value_pairs)
        
        # Update external tables
        self._insert_external_tables(extra_tables, id)

    def read_external_table(self, name, id):
        """Read row from external table."""
        con = self.connection or self._connect()
        cur = con.cursor()
        cur.execute("SELECT * FROM {} WHERE id=?".format(name), (id,))
        items = cur.fetchall()
        dictionary = dict([(item[0], item[1]) for item in items])

        if self.connection is None:
            con.close()
        return dictionary

    def _get_tables_names(self):
        """Return a list with all table names."""
        con = self.connection or self._connect()
        cur = con.cursor()
        cur.execute("SELECT name FROM sqlite_master WHERE type='table'")
        tables = cur.fetchall()
        names = [item[0] for item in tables]

        if self.connection is None:
            con.close()
        return names

    def get_external_table_names(self):
        """Return a list with all external table names."""
        from ase.db.sqlite import all_tables
        ignore = ['information', 'sqlite_sequence']
        names = self._get_tables_names()
        return list(set(names) - set(all_tables) - set(ignore))

    def _convert_tuple_to_row(self, values):
        """Convert tuples to row."""
        atoms_row = SQLite3Database._convert_tuple_to_row(self, values)

        # Now we need to update with info from the external tables
        external_tab = self.get_external_table_names()
        tables = {}
        for tab in external_tab:
            row = self.read_external_table(tab, atoms_row.id)
            tables[tab] = row
        
        if tables:
            atoms_row.__dict__.update(tables)
        return atoms_row

    def _delete(self, cur, ids, tables=None):
        from ase.db.sqlite import all_tables
        default_tables = all_tables[::-1] + self.get_external_table_names()
        tables = tables or default_tables
        SQLite3Database._delete(self, cur, ids, tables=tables)
        
