"""Use C-extensions from ase_koopmanscext.

This module defines a decorator that can be used to replace pure Python
functions with faster C-implementations from the ase_koopmans_ext module.
"""

import functools

try:
    import ase_koopmans_ext
except ImportError:
    ase_koopmans_ext = None


def cextension(func):
    if ase_koopmans_ext is None:
        return func
    cfunc = getattr(ase_koopmans_ext, func.__name__, None)
    if cfunc is None:
        return func
    functools.update_wrapper(cfunc, func)
    cfunc.__pure_python_function__ = func
    return cfunc
