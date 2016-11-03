************
Installation
************


Installing from pypi
====================

The simplest way to install ``iCount`` is from ``PyPI``. Installing is done with
one-line command::

    pip install iCount

We recommend instlling the package in `virtualenv`_. If installing the package
globally, you may need root priviliges (prefix upper command with ``sudo``)

.. _`virtualenv`:
    https://virtualenv.pypa.io/en/stable/


Installing from source
======================

If you wish to install form source, follow instructions in :doc:`Contributing
section <contributing>`.

Using docker
============

When working with docker, make sure that the docker-machine has enough memory:

    docker-machine create -d virtualbox --virtualbox-memory 65536 default
