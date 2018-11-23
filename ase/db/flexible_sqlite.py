from ase.db.sqlite import SQLite3Database
from ase.parallel import parallel_function
from ase.db.core import lock, Database

class FlexibleSqlite(SQLite3Database):
    def _create_if_not_exists(self, name, entries):
        """
        Create a new table if it does not already exists
        """
        con = self.connection or self._connect()
        types = self._get_types(entries)
        cur = con.cursor()
        cols = ",".join(k + " " + v for k, v in types.items())
        sql = "CREATE TABLE IF NOT EXISTS {} ({})".format(name, cols)
        cur.execute(sql)

        if self.connection is None:
            con.commit()
            con.close()

    def _get_types(self, entries):
        types = {}
        for k, v in entries.items():
            if isinstance(v, float):
                types[k] = "float"
            elif isinstance(v, int):
                types[k] = "integer"
            elif isinstance(v, str):
                types[k] = "text"
        return types

    def _sync_columns(self, name, entries):
        """
        Synchronize columns. (i.e. create the ones that does not exist)
        """
        con = self.connection or self._connect()
        cur = con.cursor()
        types = self._get_types(entries)

        for k, v in types.items():
            try:
                sql = "ALTER TABLE {} ADD COLUMN {} {}".format(name, k, v)
                cur.execute(sql)
            except Exception:
                # Column already exists
                pass
        
        if self.connection is None:
            con.commit()
            con.close()

    @parallel_function
    @lock
    def insert(self, name=None, entries=None, sync=True):
        """Insert a new row into a table."""
        if name is None or entries is None:
            # There is nothing to do
            return

        if sync:
            self._create_if_not_exists(name, entries)
            self._sync_columns(name, entries)
        
        con = self.connection or self._connect()
        cur = con.cursor()

        key_value = [(k, v) for k, v in entries.items()]
        if self._id_exists_in_table(cur, name, entries["id"]):
            sql = "UPDATE {} SET ".format(name)
            vals = list(item[1] for item in key_value if item[0] != "id")
            sql += ",".join(item[0]+"=?" for item in key_value if item[0] != "id")
            sql += " WHERE id=?"
            vals.append(entries["id"])
            cur.execute(sql, tuple(vals))
        else:
            sql = "INSERT INTO {} ".format(name)
            
            cols = ",".join(item[0] for item in key_value)
            vals = tuple(item[1] for item in key_value)
            sql += "({}) ".format(cols)
            sql += "VALUES (" + ",".join("?" for item in key_value) + ")"
            cur.execute(sql, vals)

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
        row = cur.fetchone()
        dictionary = None

        if row:
            column_names = [info[0] for info in cur.description]
            dictionary = dict([(name, value) for name, value in zip(column_names, row)])

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
        
