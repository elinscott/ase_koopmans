"""WSGI Flask-app for browsing a database_koopmans.

::

    +---------------------+
    | layout.html         |
    | +-----------------+ |    +--------------+
    | | search.html     | |    | layout.html  |
    | |     +           | |    | +---------+  |
    | | table.html ----------->| |row.html |  |
    | |                 | |    | +---------+  |
    | +-----------------+ |    +--------------+
    +---------------------+

You can launch Flask's local webserver like this::

    $ ase_koopmans db abc.db -w

or this::

    $ python3 -m ase_koopmans.db.app abc.db

"""

import io
import sys
from typing import Dict, Any
from pathlib import Path

from flask import Flask, render_template, request

from ase_koopmans.db import connect
from ase_koopmans.db.core import Database_koopmans
from ase_koopmans.formula import Formula
from ase_koopmans.db.web import create_key_descriptions, Session
from ase_koopmans.db.row import row2dct, AtomsRow
from ase_koopmans.db.table import all_columns


root = Path(__file__).parent.parent.parent
app = Flask(__name__, template_folder=str(root))

projects: Dict[str, Dict[str, Any]] = {}


static = root / 'ase_koopmans/db/static'
if not (static / 'jsmol/JSmol.min.js').is_file():
    print(f"""
WARNING:
    You don't have jsmol on your system.

    Download Jmol-*-binary.tar.gz from
    https://sourceforge.net/projects/jmol/files/Jmol/,
    extract jsmol.zip, unzip it and create a soft-link:

        $ tar -xf Jmol-*-binary.tar.gz
        $ unzip jmol-*/jsmol.zip
        $ ln -s $PWD/jsmol {static}/jsmol
""",
          file=sys.stderr)


@app.route('/', defaults={'project_name': 'default'})
@app.route('/<project_name>')
@app.route('/<project_name>/')
def search(project_name: str):
    """Search page.

    Contains input form for database_koopmans query and a table result rows.
    """
    session = Session(project_name)
    project = projects[project_name]
    return render_template(project['search_template'],
                           q=request.args.get('query', ''),
                           p=project,
                           session_id=session.id)


@app.route('/update/<int:sid>/<what>/<x>/')
def update(sid: int, what: str, x: str):
    """Update table of rows inside search page.

    ``what`` must be one of:

    * query: execute query in request.args (x not used)
    * limit: set number of rows to show to x
    * toggle: toggle column x
    * sort: sort after column x
    * page: show page x
    """
    session = Session.get(sid)
    project = projects[session.project_name]
    session.update(what, x, request.args, project)
    table = session.create_table(project['database_koopmans'], project['uid_key'])
    return render_template('ase_koopmans/db/templates/table.html',
                           t=table,
                           p=project,
                           s=session)


@app.route('/<project_name>/row/<uid>')
def row(project_name: str, uid: str):
    """Show details for one database_koopmans row."""
    project = projects[project_name]
    uid_key = project['uid_key']
    row = project['database_koopmans'].get('{uid_key}={uid}'
                                  .format(uid_key=uid_key, uid=uid))
    dct = project['row_to_dict_function'](row, project)
    return render_template(project['row_template'],
                           d=dct, row=row, p=project, uid=uid)


@app.route('/atoms/<project_name>/<int:id>/<type>')
def atoms(project_name: str, id: int, type: str):
    """Return atomic structure as cif, xyz or json."""
    row = projects[project_name]['database_koopmans'].get(id=id)
    a = row.toatoms()
    if type == 'cif':
        b = io.BytesIO()
        a.pbc = True
        a.write(b, 'cif', wrap=False)
        return b.getvalue(), 200, []

    fd = io.StringIO()
    if type == 'xyz':
        a.write(fd, 'xyz')
    elif type == 'json':
        con = connect(fd, type='json')
        con.write(row,
                  data=row.get('data', {}),
                  **row.get('key_value_pairs', {}))
    else:
        1 / 0

    headers = [('Content-Disposition',
                'attachment; filename="{project_name}-{id}.{type}"'
                .format(project_name=project_name, id=id, type=type))]
    txt = fd.getvalue()
    return txt, 200, headers


@app.route('/gui/<int:id>')
def gui(id: int):
    """Pop ud ase_koopmans gui window."""
    from ase_koopmans.visualize import view
    atoms = projects['default']['database_koopmans'].get_atoms(id)
    view(atoms)
    return '', 204, []


@app.route('/test')
def test():
    from pyjokes import get_joke as j
    return j()


@app.route('/robots.txt')
def robots():
    return ('User-agent: *\n'
            'Disallow: /\n'
            '\n'
            'User-agent: Baiduspider\n'
            'Disallow: /\n'
            '\n'
            'User-agent: SiteCheck-sitecrawl by Siteimprove.com\n'
            'Disallow: /\n',
            200)


def handle_query(args) -> str:
    """Converts request args to ase_koopmans.db query string."""
    return args['query']


def row_to_dict(row: AtomsRow, project: Dict[str, Any]) -> Dict[str, Any]:
    """Convert row to dict for use in html template."""
    dct = row2dct(row, project['key_descriptions'])
    dct['formula'] = Formula(Formula(row.formula).format('abc')).format('html')
    return dct


def add_project(db: Database_koopmans) -> None:
    """Add database_koopmans to projects with name 'default'."""
    all_keys = set()
    for row in db.select(columns=['key_value_pairs'], include_data=False):
        all_keys.update(row._keys)
    kd = {key: (key, '', '') for key in all_keys}
    projects['default'] = {
        'name': 'default',
        'uid_key': 'id',
        'key_descriptions': create_key_descriptions(kd),
        'database_koopmans': db,
        'row_to_dict_function': row_to_dict,
        'handle_query_function': handle_query,
        'default_columns': all_columns[:],
        'search_template': 'ase_koopmans/db/templates/search.html',
        'row_template': 'ase_koopmans/db/templates/row.html'}


if __name__ == '__main__':
    db = connect(sys.argv[1])
    add_project(db)
    app.run(host='0.0.0.0', debug=True)
