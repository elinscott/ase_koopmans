import pytest
import os
from ase_koopmans.db import connect

cmd = """
ase_koopmans -T build H | ase_koopmans -T run emt -o testase_koopmans.json &&
ase_koopmans -T build H2O | ase_koopmans -T run emt -o testase_koopmans.json &&
ase_koopmans -T build O2 | ase_koopmans -T run emt -o testase_koopmans.json &&
ase_koopmans -T build H2 | ase_koopmans -T run emt -f 0.02 -o testase_koopmans.json &&
ase_koopmans -T build O2 | ase_koopmans -T run emt -f 0.02 -o testase_koopmans.json &&
ase_koopmans -T build -x fcc Cu | ase_koopmans -T run emt -E 5,1 -o testase_koopmans.json &&
ase_koopmans -T db -v testase_koopmans.json natoms=1,Cu=1 --delete --yes &&
ase_koopmans -T db -v testase_koopmans.json "H>0" -k hydro=1,abc=42,foo=bar &&
ase_koopmans -T db -v testase_koopmans.json "H>0" --delete-keys foo"""


names = ['testase_koopmans.json',
         'testase_koopmans.db',
         'postgresql',
         'mysql',
         'mariadb']


@pytest.mark.slow
@pytest.mark.parametrize('name', names)
def test_db(name, cli):
    def count(n, *args, **kwargs):
        m = len(list(con.select(columns=['id'], *args, **kwargs)))
        assert m == n, (m, n)

    if name == 'postgresql':
        pytest.importorskip('psycopg2')
        if os.environ.get('POSTGRES_DB'):  # gitlab-ci
            name = 'postgresql://ase_koopmans:ase_koopmans@postgres:5432/testase_koopmans'
        else:
            name = os.environ.get('ASE_TEST_POSTGRES_URL')
            if name is None:
                return
    elif name == 'mysql':
        pytest.importorskip('pymysql')
        if os.environ.get('CI_PROJECT_DIR'):  # gitlab-ci
            name = 'mysql://root:ase_koopmans@mysql:3306/testase_koopmans_mysql'
        else:
            name = os.environ.get('MYSQL_DB_URL')

        if name is None:
            return
    elif name == 'mariadb':
        pytest.importorskip('pymysql')
        if os.environ.get('CI_PROJECT_DIR'):  # gitlab-ci
            name = 'mariadb://root:ase_koopmans@mariadb:3306/testase_koopmans_mysql'
        else:
            name = os.environ.get('MYSQL_DB_URL')

        if name is None:
            return

    con = connect(name)
    if 'postgres' in name or 'mysql' in name or 'mariadb' in name:
        con.delete([row.id for row in con.select()])

    cli.shell(cmd.replace('testase_koopmans.json', name))
    assert con.get_atoms(H=1)[0].magmom == 1
    count(5)
    count(3, 'hydro')
    count(0, 'foo')
    count(3, abc=42)
    count(3, 'abc')
    count(0, 'abc,foo')
    count(3, 'abc,hydro')
    count(0, foo='bar')
    count(1, formula='H2')
    count(1, formula='H2O')
    count(3, 'fmax<0.1')
    count(1, '0.5<mass<1.5')
    count(5, 'energy')

    id = con.reserve(abc=7)
    assert con[id].abc == 7

    for key in ['calculator', 'energy', 'abc', 'name', 'fmax']:
        count(6, sort=key)
        count(6, sort='-' + key)

    cli.shell('ase_koopmans -T gui --terminal -n 3 {}'.format(name))

    con.delete([id])
