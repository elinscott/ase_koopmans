import re
from typing import List, Tuple, Dict

from flask import flash

from ase.db.core import default_key_descriptions, Database
from ase.db.table import Table, all_columns


class Session:
    next_id = 1
    sessions: Dict[int, 'Session'] = {}

    def __init__(self):
        self.id = Session.next_id
        Session.next_id += 1

        Session.sessions[self.id] = self
        if len(Session.sessions) > 1000:
            # Forget old sessions:
            for id in sorted(Session.sessions)[:200]:
                del Session.sessions[id]

        self.columns = None
        self.nrows = None
        self.page = 0
        self.limit = 25
        self.sort = ''
        self.query = ''
        self.project = 'default'
        self.create_table_functions = None
        self.unique_key = 'id'

    def __str__(self):
        return str(self.__dict__)

    @staticmethod
    def get(id: int) -> 'Session':
        if id in Session.sessions:
            return Session.sessions[id]
        return Session()

    def update(self,
               what: str,
               x: str,
               query: str,
               default_columns: List[str]) -> None:

        if self.columns is None:
            self.columns = default_columns[:]

        if what == 'query':
            self.query = query
            self.nrows = None

        elif what == 'sort':
            if x == self.sort:
                self.sort = '-' + x
            elif '-' + x == self.sort:
                self.sort = 'id'
            else:
                self.sort = x
            self.page = 0

        elif what == 'limit':
            self.limit = int(x)
            self.page = 0

        elif what == 'page':
            self.page = int(x)

        elif what == 'toggle':
            column = x
            if column == 'reset':
                self.columns = default_columns[:]
            else:
                if column in self.columns:
                    self.columns.remove(column)
                    if column == self.sort.lstrip('-'):
                        self.sort = 'id'
                        self.page = 0
                else:
                    self.columns.append(column)

    @property
    def row1(self) -> int:
        return self.page * self.limit + 1

    @property
    def row2(self) -> int:
        return min((self.page + 1) * self.limit, self.nrows)

    def paginate(self) -> List[Tuple[int, str]]:
        """Helper function for pagination stuff."""
        npages = (self.nrows + self.limit - 1) // self.limit
        p1 = min(5, npages)
        p2 = max(self.page - 4, p1)
        p3 = min(self.page + 5, npages)
        p4 = max(npages - 4, p3)
        pgs = list(range(p1))
        if p1 < p2:
            pgs.append(-1)
        pgs += list(range(p2, p3))
        if p3 < p4:
            pgs.append(-1)
        pgs += list(range(p4, npages))
        pages = [(self.page - 1, 'previous')]
        for p in pgs:
            if p == -1:
                pages.append((-1, '...'))
            elif p == self.page:
                pages.append((-1, str(p + 1)))
            else:
                pages.append((p, str(p + 1)))
        nxt = min(self.page + 1, npages - 1)
        if nxt == self.page:
            nxt = -1
        pages.append((nxt, 'next'))
        return pages

    def create_table(self,
                     db: Database) -> Table:

        if self.create_table_function is not None:
            return self.create_table_function(self, db, self.unique_key)

        query = self.query
        if self.nrows is None:
            try:
                self.nrows = db.count(query)
            except (ValueError, KeyError) as e:
                error = ', '.join(['Bad query'] + list(e.args))
                flash(error)
                query = 'id=0'  # this will return no rows
                self.nrows = 0

        table = Table(db, self.unique_key)
        table.select(query, self.columns, self.sort,
                     self.limit, offset=self.page * self.limit)
        table.format()
        table.addcolumns = sorted(column for column in
                                  all_columns + table.keys
                                  if column not in table.columns)
        return table


KeyDescriptions = Dict[str, Tuple[str, str, str]]  # type-hint shortcut


def create_key_descriptions(kd: KeyDescriptions) -> KeyDescriptions:
    kd = kd.copy()
    kd.update(default_key_descriptions)

    # Long description may be missing:
    for key, (short, long, unit) in kd.items():
        if not long:
            kd[key] = (short, short, unit)

    sub = re.compile(r'`(.)_(.)`')
    sup = re.compile(r'`(.*)\^\{?(.*?)\}?`')

    # Convert LaTeX to HTML:
    for key, value in kd.items():
        short, long, unit = value
        unit = sub.sub(r'\1<sub>\2</sub>', unit)
        unit = sup.sub(r'\1<sup>\2</sup>', unit)
        unit = unit.replace(r'\text{', '').replace('}', '')
        kd[key] = (short, long, unit)

    return kd
