import os
from ase.db import connect
from ase import Atoms
from ase.test import must_raise

DB_NAMES = ["test_ext_tables.db", "postgresql"]

def get_db_name(name):
    if name == 'postgresql':
        if os.environ.get('POSTGRES_DB'):  # gitlab-ci
            name = 'postgresql://ase:ase@postgres:5432/testase'
        else:
            name = os.environ.get('ASE_TEST_POSTGRES_URL')
    return name

def test_create_and_delete_ext_tab(db_name):
    ext_tab = ["tab1", "tab2", "tab3"]
    atoms = Atoms()
    db = connect(name)
    db.write(atoms)

    for tab in ext_tab:
        db.create_table_if_not_exists(tab, "INTEGER")
    assert sorted(db.get_external_table_names()) == ext_tab

    db.delete_external_table("tab1")
    assert sorted(db.get_external_table_names()) == ["tab2", "tab3"]

def test_insert_in_external_tables(db_name):
    atoms = Atoms()
    db = connect(db_name)

    # Now a table called insert_tab with schema datatype REAL should
    # be created
    uid = db.write(atoms, external_tables={"insert_tab": {"rate": 1.0, "rate1": -2.0}})

    db.delete([uid])

    # Hack: retrieve the connection
    con = db._connect()
    cur = con.cursor()

    sql = "SELECT * FROM insert_tab WHERE ID=?"
    cur.execute(sql, (uid,))

    entries = [x for x in cur.fetchall()]
    if db.connection is None:
        con.close()
    assert not entries

    # Make sure that there are now entries in the 
    # external table with current uid

    # Try to insert something that should not pass
    # i.e. string value into the same table
    with must_raise(ValueError):
        db.write(atoms, external_tables={"insert_tab": {"rate": "something"}})

def test_extract_from_table(db_name):
    atoms = Atoms()
    db = connect(db_name)
    uid = db.write(atoms, external_tables={"insert_tab": {"rate": 12.0, "rate1": -10.0}})

    row = db.get(id=uid)
    assert abs(row["insert_tab"]["rate"] - 12.0) < 1E-8
    assert abs(row["insert_tab"]["rate1"] + 10.0) < 1E-8

def test_write_atoms_row(db_name):
    atoms = Atoms()
    db = connect(db_name)
    uid = db.write(atoms, external_tables={"insert_tab": {"rate": 12.0, "rate1": -10.0}, 
                                           "another_tab": {"somevalue": 1.0}})
    row = db.get(id=uid)

    # Hack: Just change the unique ID 
    row["unique_id"] = "uniqueIDTest"
    db.write(row)



for db_name in DB_NAMES:
    name = get_db_name(db_name)

    if name is None:
        continue
    test_create_and_delete_ext_tab(db_name)
    test_insert_in_external_tables(db_name)
    test_extract_from_table(db_name)
    test_write_atoms_row(db_name)
    os.remove(db_name)

