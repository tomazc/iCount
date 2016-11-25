************
Installation
************


Installing from pypi
====================

The simplest way to install **iCount** is from `PyPI`_. Installing is done with
one-line command::

    pip install iCount

We recommend installing the package in `virtualenv`_. If installing the package
globally, you may need root priviliges (prefix upper command with ``sudo``)

.. _`virtualenv`:
    https://virtualenv.pypa.io/en/stable/

.. _`PyPI`:
    https://pypi.python.org/pypi


Installing from source
======================

If you wish to install form source, follow instructions in section
:doc:`Contributing <contributing>`.


Docker
======

When working with docker, make sure that the docker-machine has enough memory to run STAR and
associated programs, *e.g.*, at least 64 GB::

    docker-machine create -d virtualbox --virtualbox-memory 65536 default
