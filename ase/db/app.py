"""WSGI Flask-app for browsing a database.

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

    $ ase db abc.db -w

or this::

    $ python3 -m ase.db.app abc.db

"""

import io
import sys
from typing import Dict, Any

from flask import Flask, render_template, request

from ase.db import connect
from ase.formula import Formula
from ase.db.web import create_key_descriptions, Session
from ase.db.row import row2dct
from ase.db.table import all_columns


app = Flask(__name__)

projects = {}  # type: Dict[str, Dict[str, Any]]


@app.route('/', defaults={'project_name': 'default'})
@app.route('/<project_name>/')
@app.route('/<project_name>')
def search(project_name: str):
    session = Session(project_name)
    project = projects[project_name]
    return render_template(project['search_template'],
                           p=project,
                           session_id=session.id)


@app.route('/update/<int:sid>/<what>/<x>/')
def update(sid, what, x):
    session = Session.get(sid)
    project = projects[session.project_name]
    session.update(what, x, request.args, project)
    table = session.create_table(project['database'], project['uid_key'])
    return render_template('table.html',
                           t=table,
                           p=project,
                           s=session)


@app.route('/<project_name>/row/<uid>')
def row(project_name: str, uid: str):
    project = projects[project_name]
    uid_key = project['uid_key']
    row = project['database'].get('{uid_key}={uid}'
                                  .format(uid_key=uid_key, uid=uid))
    dct = project['row_to_dict_function'](row, project)
    return render_template(project['row_template'],
                           d=dct, row=row, p=project)


@app.route('/atoms/<project_name>/<int:id>/<type>')
def atoms(project_name, id, type):
    row = projects[project_name]['database'].get(id=id)
    a = row.toatoms()
    if type == 'cif':
        # fd = io.BytesIO()
        # a.write(fd, 'cif')
        # return fd.getvalue(), 200, []

        a.write('x.cif')
        from flask import send_from_directory
        return send_from_directory('.', 'x.cif')

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
    from ase.visualize import view
    atoms = projects['default']['database'].get_atoms(id)
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


def handle_query(args):
    return args['query']


def row_to_dict(row, project):
    dct = row2dct(row, project['key_descriptions'])
    dct['formula'] = Formula(Formula(row.formula).format('abc')).format('html')
    return dct


def init(db):
    all_keys = set()
    for row in db.select(columns=['key_value_pairs'], include_data=False):
        all_keys.update(row._keys)
    kd = {key: (key, '', '') for key in all_keys}
    projects['default'] = {
        'name': 'default',
        'uid_key': 'id',
        'key_descriptions': create_key_descriptions(kd),
        'database': db,
        'row_to_dict_function': row_to_dict,
        'handle_query_function': handle_query,
        'default_columns': all_columns[:],
        'search_template': 'search.html',
        'row_template': 'row.html'}


if __name__ == '__main__':
    db = connect(sys.argv[1])
    init(db)
    app.run(host='0.0.0.0', debug=True)
