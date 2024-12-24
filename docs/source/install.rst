Install
=======

ModiFinder requires Python 3.11 or above. 

This turorial assumes you have already setup the Python environment you want to work in and intend to install ``modifinder`` inside of it.  If you want
to create and work with Python virtual environments, please follow instructions
on `venv <https://docs.python.org/3/library/venv.html>`_ and `virtual
environments <http://docs.python-guide.org/en/latest/dev/virtualenvs/>`_ or use `Conda/Mamba <https://github.com/mamba-org/mamba>`_ to create a new environment.

First, make sure you have the latest version of ``pip`` installed. `Pip documentation
<https://pip.pypa.io/en/stable/installing/>`_ 

Install the released version
-----------------------------
Install the current release of ``modifinder`` with ``pip``::

    $ pip install modifinder

To upgrade to a newer release use the ``--upgrade`` flag::

    $ pip install --upgrade modifinder

Alternatively, you can manually download ``modifinder`` from
`GitHub <https://github.com/Wang-Bioinformatics-Lab/ModiFinder_base>`_
and install it from the source code.  To do this, unpack the source code and
run the following from the top-level source directory using the Terminal::

    $ pip install .