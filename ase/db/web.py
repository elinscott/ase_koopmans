import re
from typing import List, Tuple, Dict

from ase.db.core import default_key_descriptions


def paginate(page: int, nrows: int, limit: int) -> List[Tuple[int, str]]:
    """Helper function for pagination stuff."""
    npages = (nrows + limit - 1) // limit
    p1 = min(5, npages)
    p2 = max(page - 4, p1)
    p3 = min(page + 5, npages)
    p4 = max(npages - 4, p3)
    pgs = list(range(p1))
    if p1 < p2:
        pgs.append(-1)
    pgs += list(range(p2, p3))
    if p3 < p4:
        pgs.append(-1)
    pgs += list(range(p4, npages))
    pages = [(page - 1, 'previous')]
    for p in pgs:
        if p == -1:
            pages.append((-1, '...'))
        elif p == page:
            pages.append((-1, str(p + 1)))
        else:
            pages.append((p, str(p + 1)))
    nxt = min(page + 1, npages - 1)
    if nxt == page:
        nxt = -1
    pages.append((nxt, 'next'))
    return pages


def create_key_descriptions(db) -> Dict[str, Tuple[str, str, str]]:
    kd = default_key_descriptions.copy()

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

    all_keys = set()
    for row in db.select(columns=['key_value_pairs'], include_data=False):
        all_keys.update(row._keys)
    for key in all_keys:
        kd[key] = (key, key, '')

    return kd
