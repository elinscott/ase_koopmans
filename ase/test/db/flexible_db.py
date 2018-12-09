from ase.db import connect
from ase.build import bulk
import os

def read_write_no_tables():
    """Test that we can do things as before."""
    db_name = "read_write_no_tables.db"
    db = connect(db_name, type="flexsqlite")

    atoms = bulk("Al")*(4, 4, 4)
    db.write(atoms, tag="pure_al")

    atoms_si = atoms.copy()
    atoms_si[0].symbol = "Si"
    db.write(atoms_si, tag="pure_si")

    atoms_si_db = db.get(tag="pure_si").toatoms()
    assert atoms_si_db == atoms_si

    atoms_al_db = db.get(tag="pure_al").toatoms()
    assert atoms_al_db == atoms
    os.remove(db_name)

def read_write_table():
    """Test if we can write an atoms object with an attached table."""

    db_name = "read_write_table_test.db"
    db = connect(db_name, type="flexsqlite")
    atoms = bulk("Al")*(4, 4, 4)
    thermo = {"temperature": 100, "Cv": 0.05, "internal_enregy": -20.0, "chemical_potential": -0.05}

    db.write(atoms, tables={"thermo": thermo})

    # Verify that it works to add a new entry with a different set of 
    # columns
    thermo["free_energy"] = 2.0
    db.write(atoms, tables={"thermo": thermo})

    db.write(atoms, tables={"test_table": {"some_col": 1.0}})
    external = db.get_external_table_names()
    assert external == set(["test_table", "thermo"])

    row = db.get(id=2)
    for k, v in thermo.items():
        if row["thermo"][k] is not None and k != "id":
            assert abs(row["thermo"][k] - v) < 1E-8

    # Try to update the table
    thermo2 = {"temperature": 200, "Cv": 0.035, "internal_enregy": -25.0, "chemical_potential": -0.1}
    db.update(1, tables={"thermo": thermo2})

    # Read back and verify
    row = db.get(id=1)
    for k, v in thermo2.items():
        if row["thermo"][k] is not None and k != "id":
            assert abs(row["thermo"][k] - v) < 1E-8

    # Try do delete
    del db[1]

    # Verify that there are no entries with ID 1 in the external tables
    import sqlite3
    con = sqlite3.connect(db_name)
    cur = con.cursor()
    cur.execute("SELECT id FROM thermo WHERE id=1")
    exists_in_thermo = cur.fetchone() is not None
    con.close()
    assert not exists_in_thermo

    # Try to insert many entries within a context manager
    with db:
        for _ in range(10):
            db.write(atoms, tables={"thermo1": thermo, "thermo2": thermo2})
    os.remove(db_name)


def write_atoms_row():
    db_name = "write_atoms_row.db"
    db = connect(db_name, type="flexsqlite")
    
    atoms = bulk("Ta")*(4, 4, 4)
    uid = db.write(atoms, tables={"some_table": {"val1": 0.0, "val2": 1.0}})

    # Retrieve the row from the database
    row = db.get(id=uid)

    # Hack: Update the __dict__ entry of the row
    row.__dict__["some_table"]["val2"] = 6.0
    row.__dict__["unique_id"] = "1e45a76"


    # Update the database
    uid = db.write(row)

    # Again get the row and verify that the table
    # has to updated values
    row2 = db.get(id=uid)
    assert abs(row2["some_table"]["val1"]) < 1E-6
    assert abs(row2["some_table"]["val2"] - 6.0) < 1E-6
    os.remove(db_name)


read_write_no_tables()
read_write_table()
write_atoms_row()




    