.. _download_and_install:

============
Installation
============

Requirements
============

* Python_ 3.6 or newer
* NumPy_ 1.11 or newer (base N-dimensional array package)
* SciPy_ 0.18 or newer (library for scientific computing)

Optional but strongly recommended:

* Matplotlib_ 2.0.0 or newer for plotting
* :mod:`tkinter` for :mod:`ase.gui`

Optional:

* Flask_ for :mod:`ase.db` web-interface
* pytest_ 3.6.1 or newer for running tests
* pytest-xdist_ 1.22.1 or newer for running tests in parallel
* spglib_ for certain symmetry-related features

.. _Python: https://www.python.org/
.. _NumPy: https://docs.scipy.org/doc/numpy/reference/
.. _SciPy: https://docs.scipy.org/doc/scipy/reference/
.. _Matplotlib: https://matplotlib.org/
.. _Flask: https://palletsprojects.com/p/flask/
.. _PyPI: https://pypi.org/project/ase
.. _PIP: https://pip.pypa.io/en/stable/
.. _pytest: https://pypi.org/project/pytest/
.. _pytest-xdist: https://pypi.org/project/pytest-xdist/
.. _spglib: https://pypi.org/project/spglib/

Installation using system package managers
==========================================

Linux
-----

Major GNU/Linux distributions (including Debian and Ubuntu derivatives,
Arch, Fedora, Red Hat and CentOS) have a ``python-ase`` package
available that you can install on your system. This will manage
dependencies and make ASE available for all users.

.. note::
   Depending on the distribution, this may not be the latest
   release of ASE.

Max OSX (Homebrew)
------------------

The old version of Python included in Mac OSX is incompatible with ASE
and does not include the pip_ package manager.

Before installing ASE with ``pip`` as described in the next section, Mac
users need to install an appropriate Python version.  One option is
to use the Homebrew_ package manager, which provides an up-to-date version
of Python 3 including ``pip`` and the tkinter graphical interface bindings::

  $ brew install python

For more information about the quirks of brewed Python see this guide_.

.. _Homebrew: http://brew.sh

.. _guide: https://docs.brew.sh/Homebrew-and-Python


.. index:: pip
.. _pip installation:


Installation using pip
======================

.. highlight:: bash

The simplest way to install ASE is to use pip_ which will automatically get
the source code from PyPI_::

    $ pip install --upgrade --user ase

If you intend to run the tests, use::

    $ pip install --upgrade --user ase[test]

This will install ASE in a local folder where Python can
automatically find it (``~/.local`` on Unix, see here_ for details).  Some
:ref:`cli` will be installed in the following location:

=================  ============================
Unix and Mac OS X  ``~/.local/bin``
Homebrew           ``~/Library/Python/X.Y/bin``
Windows            ``%APPDATA%/Python/Scripts``
=================  ============================

Make sure you have that path in your :envvar:`PATH` environment variable.

Now you should be ready to use ASE, but before you start, you may
wish to `run the tests`_ as described below.


.. note::

    If your OS doesn't have ``numpy``, ``scipy`` and ``matplotlib`` packages
    installed, you can install them with::

        $ pip install --upgrade --user numpy scipy matplotlib


.. _here: https://docs.python.org/3/library/site.html#site.USER_BASE


.. _download:

Installation from source
========================

As an alternative to ``pip``, you can also get the source from a tar-file or
from Git.

:Tar-file:

    You can get the source as a `tar-file <http://xkcd.com/1168/>`__ for the
    latest stable release (ase-3.19.1.tar.gz_) or the latest
    development snapshot (`<snapshot.tar.gz>`_).

    Unpack and make a soft link::

        $ tar -xf ase-3.19.1.tar.gz
        $ ln -s ase-3.19.1 ase

    Here is a `list of tarballs <https://pypi.org/simple/ase/>`__.

:Git clone:

    Alternatively, you can get the source for the latest stable release from
    https://gitlab.com/ase/ase like this::

        $ git clone -b 3.19.1 https://gitlab.com/ase/ase.git

    or if you want the development version::

        $ git clone https://gitlab.com/ase/ase.git

:Pip:

    install git master directly with pip::

        $ pip install --upgrade git+https://gitlab.com/ase/ase.git@master

    The ``--upgrade`` ensures that you always reinstall even if the version
    number hasn't changed.


Add ``~/ase`` to your :envvar:`PYTHONPATH` environment variable and add
``~/ase/bin`` to :envvar:`PATH` (assuming ``~/ase`` is where your ASE
folder is).  Alternatively, you can install the code with ``python setup.py
install --user`` and add ``~/.local/bin`` to the front of your :envvar:`PATH`
environment variable (if you don't already have that).

Finally, please `run the tests`_.

.. note::

    We also have Git-tags for older stable versions of ASE.
    See the :ref:`releasenotes` for which tags are available.  Also the
    dates of older releases can be found there.


.. _ase-3.19.1.tar.gz: https://pypi.org/packages/source/a/ase/ase-3.19.1.tar.gz


Environment variables
=====================

.. envvar:: PATH

    Colon-separated paths where programs can be found.

.. envvar:: PYTHONPATH

    Colon-separated paths where Python modules can be found.

Set these permanently in your :file:`~/.bashrc` file::

    $ export PYTHONPATH=<path-to-ase-package>:$PYTHONPATH
    $ export PATH=<path-to-ase-command-line-tools>:$PATH

or your :file:`~/.cshrc` file::

    $ setenv PYTHONPATH <path-to-ase-package>:${PYTHONPATH}
    $ setenv PATH <path-to-ase-command-line-tools>:${PATH}

.. note::

   If running on Mac OSX: be aware that terminal sessions will
   source :file:`~/.bash_profile` by default and not
   :file:`~/.bashrc`. Either put any ``export`` commands into
   :file:`~/.bash_profile` or source :file:`~/.bashrc` in all Bash
   sessions by adding

   ::

      if [ -f ${HOME}/.bashrc ]; then
      source ${HOME}/.bashrc
      fi

   to your :file:`~/.bash_profile`.


.. index:: test
.. _running tests:
.. _run the tests:

Test your installation
======================

Before running the tests, make sure you have set your :envvar:`PATH`
environment variable correctly as described in the relevant section above.
Run the tests like this::

    $ ase test  # takes 1 min.

and send us the output if there are failing tests.
