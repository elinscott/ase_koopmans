"""WSGI Flask-app for browsing a database.

You can launch Flask's local webserver like this::

    $ ase db abc.db -w

For a real webserver, you need to set the $ASE_DB_APP_CONFIG environment
variable to point to a configuration file like this::

    ASE_DB_NAMES = ['/path/to/db-file/project1.db',
                    'postgresql://user:pw@localhost:5432/project2']
    ASE_DB_TEMPLATES = '...'
    ASE_DB_TMPDIR = '...'
    ASE_DB_DOWNLOAD = False  # or True

Start with something like::

    twistd web --wsgi=ase.db.app.app --port=8000

"""

import io
import sys

from flask import Flask, render_template, request, flash

import ase.db
from ase.db.web import (create_key_descriptions, create_table, Session,
                        row2things)
from ase.db.row import AtomsRow
from ase.db.table import all_columns


app = Flask(__name__)

databases = {}
key_descriptions = {}


def get_row(pid: str) -> AtomsRow:
    project, _, s = pid.rpartition('-')
    id = int(s)
    return databases[project][id]


@app.route('/')
def index():
    session_id = int(request.args.get('x', '0'))
    session = Session.get(session_id)
    session.update(request.args.get('page', None),
                   request.args.get('limit', 0),
                   request.args.get('sort', ''),
                   request.args.get('toggle', ''),
                   all_columns)
    db = databases['']
    table, error = create_table(db, session)
    if error:
        flash(error)
    return render_template('table.html',
                           t=table,
                           kd=key_descriptions,
                           s=session)


@app.route('/row/<int:id>')
def row(id):
    db = databases['']
    row = db.get(id=id)
    things = row2things(row, key_descriptions)
    return render_template('summary.html', t=things)


@app.route('/atoms/<pid>/<type>')
def atoms(pid, type):
    row = get_row(pid)
    a = row.toatoms()
    if type == 'cif':
        fd = io.StringIO()
        a.write(fd, 'cif')
        return fd.getvalue(), 200, []

    fd = io.StringIO()
    if type == 'xyz':
        a.write(fd, 'cif')
    elif type == 'json':
        con = ase.db.connect(fd, type='json')
        con.write(row,
                  data=row.get('data', {}),
                  **row.get('key_value_pairs', {}))
    else:
        1 / 0

    headers = [('Content-Disposition',
                f'attachment; filename="{pid}.{type}"')]
    txt = fd.getvalue()
    return txt, 200, headers


@app.route('/gui/<int:id>')
def gui(id: int):
    from ase.visualize import view
    db = databases['']
    atoms = db.get_atoms(id)
    view(atoms)
    return '', 204, []


@app.route('/test')
def test():
    import pyjokes as j
    return j.get_joke('en')


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


def main(db):
    databases[''] = db
    key_descriptions.update(create_key_descriptions(db))


if __name__ == '__main__':
    db = ase.db.connect(sys.argv[1])
    main(db)
    app.run(host='0.0.0.0', debug=True)
