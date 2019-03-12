import sys
import numpy as np

py3 = sys.version_info[0] == 3

inplace_methods = ['__iadd__', '__imul__', '__ipow__', '__isub__',
                   '__itruediv__', '__imatmul__']

forward_methods = ['__abs__', '__add__', '__contains__', '__eq__',
                   '__ge__', '__getitem__', '__gt__', '__hash__',
                   '__iter__', '__le__', '__len__', '__lt__',
                   '__mul__', '__ne__', '__neg__', '__pos__',
                   '__pow__', '__radd__', '__rmul__', '__rpow__',
                   '__rsub__', '__rtruediv__', '__setitem__',
                   '__sub__', '__truediv__']

if hasattr(np.ndarray, '__matmul__'):
    forward_methods += ['__matmul__', '__rmatmul__']

forward_methods += ['all', 'any', 'diagonal', 'dot', 'mean', 'sum',
                    'ravel', 'tolist',
                    'transpose', 'tofile', 'tobytes', 'tostring']


def forward_inplace_call(name):
    arraymeth = getattr(np.ndarray, name)
    def f(self, obj):
        a = self.__array__()
        arraymeth(a, obj)
        return self
    # use update_wrapper()?
    f.__name__ = name
    f.__qualname__ = name
    return f


def forward_call(name):
    arraymeth = getattr(np.ndarray, name)
    def f(self, *args, **kwargs):
        a = self.__array__()
        return arraymeth(a, *args, **kwargs)
    f.__name__ = name
    f.__qualname__ = name
    return f


def arraylike(cls):
    for name in inplace_methods:
        if hasattr(np.ndarray, name) and not hasattr(cls, name):
            meth = forward_inplace_call(name)
        setattr(cls, name, meth)
    for name in forward_methods:
        assert hasattr(np.ndarray, name), name
        if hasattr(np.ndarray, name):  # and not hasattr(cls, name):
            meth = forward_call(name)
            setattr(cls, name, meth)
    return cls
