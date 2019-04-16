"""Use C-extensions from asecext.

This module defines a decorator that can be used to replace pure Python
functions with faster C-implementations from the ase_ext module.
"""

try:
    from ase_ext import cextension
except ImportError:
    # No ase_ext module
    def cextension(func):
        """Just return the pure Python function."""
        return func
