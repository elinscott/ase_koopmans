"""Use C-extensions from asecext.

This module defines a decorator that can be used to replace pure Python
functions with faster C-implementations from the asecext module.
"""

try:
    from asecext import cextension
except ImportError:
    # No asecext module
    def cextension(func):
        """Just return the pure Python function."""
        return func
