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

import collections
import functools
import io
import os
import os.path as op
import sys
import tempfile

from flask import Flask, render_template, request, send_from_directory, flash

import ase.db
from ase.db.web import paginate, create_key_descriptions
from ase.db.plot import atoms2png
from ase.db.row import AtomsRow
from ase.db.summary import Summary
from ase.db.table import Table, all_columns
from ase.visualize import view
from ase import Atoms
from ase.calculators.calculator import kptdensity2monkhorstpack


# Every client-connetions gets one of these tuples:
Connection = collections.namedtuple(
    'Connection',
    ['query',  # query string
     'nrows',  # number of rows matched
     'page',  # page number
     'columns',  # what columns to show
     'sort',  # what column to sort after
     'limit'])  # number of rows per page

app = Flask(__name__)

next_con_id = 1
connections = {}
databases = {}
key_descriptions = {}

tmpdir = tempfile.mkdtemp(prefix='ase-db-app-')  # used to cache png-files


def get_row(pid: str) -> AtomsRow:
    project, _, s = pid.rpartition('-')
    id = int(s)
    return databases[project][id]


def create_table(args, db, unique_key='id'):
    global next_con_id

    con_id = int(args.get('x', '0'))

    if con_id in connections:
        query, nrows, page, columns, sort, limit = connections[con_id]

    if con_id not in connections:
        # Give this connetion a new id:
        con_id = next_con_id
        next_con_id += 1
        query = ['', {}, '']
        nrows = None
        page = 0
        columns = None
        sort = 'id'
        limit = 25

    if columns is None:
        columns = list(all_columns)

    if 'sort' in args:
        column = args['sort']
        if column == sort:
            sort = '-' + column
        elif '-' + column == sort:
            sort = 'id'
        else:
            sort = column
        page = 0
    elif 'limit' in args:
        limit = int(args['limit'])
        page = 0
    elif 'page' in args:
        page = int(args['page'])

    if 'toggle' in args:
        column = args['toggle']
        if column == 'reset':
            columns = list(all_columns)
        else:
            if column in columns:
                columns.remove(column)
                if column == sort.lstrip('-'):
                    sort = 'id'
                    page = 0
            else:
                columns.append(column)

    okquery = query

    if nrows is None:
        try:
            print(query, db)
            nrows = db.count(query[2])
        except (ValueError, KeyError) as e:
            flash(', '.join(['Bad query'] + list(e.args)))
            okquery = ('', {}, 'id=0')  # this will return no rows
            nrows = 0

    table = Table(db, unique_key)
    table.select(okquery[2], columns, sort, limit, offset=page * limit)

    con = Connection(query, nrows, page, columns, sort, limit)
    connections[con_id] = con

    if len(connections) > 1000:
        # Forget old connections:
        for cid in sorted(connections)[:200]:
            del connections[cid]

    table.format()

    addcolumns = sorted(column for column in all_columns + table.keys
                        if column not in table.columns)

    return (table, con_id, addcolumns, paginate(page, nrows, limit))


@app.route('/')
def index():
    db = databases['']
    table, con_id, addcolumns, pages = create_table(request.args, db)
    con = connections[con_id]
    print(table.columns)
    print(key_descriptions)
    return render_template('table.html',
                           t=table,
                           kd=key_descriptions,
                           con=con,
                           x=con_id,
                           addcolumns=addcolumns,
                           pages=pages,
                           row1=con.page * con.limit + 1,
                           row2=min((con.page + 1) * con.limit, con.nrows))


@app.route('/image/<pid>')
def image(pid):
    path = tmpdir / pid + '.png'
    if not path.is_file():
        atoms = get_row(pid).toatoms()
        atoms2png(atoms, path)

    return send_from_directory(tmpdir, path.name)


@app.route('/<project>/cif/<name>')
def cif(project, name):
    id = int(name[:-4])
    name = project + '-' + name
    path = op.join(tmpdir, name)
    if not op.isfile(path):
        db = databases[project]
        atoms = db.get_atoms(id)
        atoms.write(path)
    return send_from_directory(tmpdir, name)


@app.route('/<project>/plot/<uid>/<png>')
def plot(project, uid, png):
    png = project + '-' + uid + '-' + png
    return send_from_directory(tmpdir, png)


@app.route('/<project>/gui/<int:id>')
def gui(project, id):
    db = databases[project]
    atoms = db.get_atoms(id)
    view(atoms)
    return '', 204, []


@app.route('/row/<int:id>')
def row(id):
    db = databases['']
    row = db.get(id=id)
    s = Summary(row)
    atoms = Atoms(cell=row.cell, pbc=row.pbc)
    n1, n2, n3 = kptdensity2monkhorstpack(atoms,
                                          kptdensity=1.8,
                                          even=False)
    return render_template('summary.html',
                           s=s,
                           id=id,
                           n1=n1,
                           n2=n2,
                           n3=n3,
                           back=True,
                           md=db.meta)


def tofile(project, query, type, limit=0):
    fd, name = tempfile.mkstemp(suffix='.' + type)
    con = ase.db.connect(name, use_lock_file=False)
    db = databases[project]
    for row in db.select(query, limit=limit):
        con.write(row,
                  data=row.get('data', {}),
                  **row.get('key_value_pairs', {}))
    os.close(fd)
    data = open(name, 'rb').read()
    os.unlink(name)
    return data


def download(f):
    @functools.wraps(f)
    def ff(*args, **kwargs):
        text, name = f(*args, **kwargs)
        if name is None:
            return text
        headers = [('Content-Disposition',
                    'attachment; filename="{}"'.format(name)),
                   ]  # ('Content-type', 'application/sqlite3')]
        return text, 200, headers
    return ff


@app.route('/<project>/xyz/<int:id>')
@download
def xyz(project, id):
    fd = io.StringIO()
    from ase.io.xyz import write_xyz
    db = databases[project]
    write_xyz(fd, db.get_atoms(id))
    data = fd.getvalue()
    return data, '{}.xyz'.format(id)


def test():
    import pyjokes as j
    return j.get_joke('en')


@app.route('/<project>/json/<int:id>')
@download
def json1(project, id):
    if project not in databases:
        return 'No such project: ' + project, None
    data = tofile(project, id, 'json')
    return data, '{}.json'.format(id)


@app.route('/<project>/sqlite/<int:id>')
@download
def sqlite1(project, id):
    if project not in databases:
        return 'No such project: ' + project, None
    data = tofile(project, id, 'db')
    return data, '{}.db'.format(id)


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


if __name__ == '__main__':
    db = ase.db.connect(sys.argv[1])
    databases[''] = db
    key_descriptions.update(create_key_descriptions(db))
    app.run(host='0.0.0.0', debug=True)
