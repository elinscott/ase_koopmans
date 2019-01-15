import os
from ase.db import connect
from ase import Atoms

DB_NAMES = ["test_ext_tables.db", "postgresql"]

def get_db_name(name):
    if name == 'postgresql':
        if os.environ.get('POSTGRES_DB'):  # gitlab-ci
            name = 'postgresql://ase:ase@postgres:5432/testase'
        else:
            name = os.environ.get('ASE_TEST_POSTGRES_URL')
    return name

def test_create_and_delete_ext_tab():
    ext_tab = ["tab1", "tab2", "tab3"]
    atoms = Atoms()
    for db_name in DB_NAMES:
        name = get_db_name(db_name)

        if name is None:
            continue
        db = connect(name)
        db.write(atoms)

        for tab in ext_tab:
            db.create_table_if_not_exists(tab, "INTEGER")
        assert sorted(db.get_external_table_names()) == ext_tab

        db.delete_external_table("tab1")
        assert sorted(db.get_external_table_names()) == ["tab2", "tab3"]

test_create_and_delete_ext_tab()

