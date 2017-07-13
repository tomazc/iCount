************
Installation
************


Running with Docker
===================

`Docker`_ provides an easy way to run iCount in a working environment that is completely separated
from your host machine.

A ready-to-use image of iCount is available at `Docker Hub`_. After
:ref:`Installing Docker <installing-docker>` on your machine, you can run iCount by issuing the
following command::

    docker run -i -t tomazc/icount


If you want to build an image from source, then change to the source folder and issue the
following command::

    docker build -t icountsrc .


You can then run the freshly built image::

    docker run -i -t icountsrc


To make all the files and results persistent, even after the container stops, create a folder
and mount it as a `data volume`_::

    mkdir `pwd`/storage_docker
    docker run -i -t -v `pwd`/storage_docker:/home/icuser/storage icountsrc

.. note::
    Make sure to create a local folder and provide the path to it. The example above uses a path
    that may not be applicable to your computer. Both, path to the folder on the host machine and
    path within the container (``/home/icuser/storage``), must be absolute.

If you are developing iCount, then it makes sense to mount the source folder as a volume into the
docker container. Make sure to change into the source folder and then issue::

    docker run -i -t -v `pwd`:/home/icuser/iCount_src \
    -v `pwd`/storage_docker:/home/icuser/storage icountsrc

This setup will behave similarly to the :doc:`Installing from source <contributing>` setup.
All changes to the source made on your host computer, will be immediately available in the
running container.


.. _`Docker`:
    https://www.docker.com

.. _`Docker Hub`:
    https://hub.docker.com/r/tomazc/icount/

.. _`data volume`:
    https://docs.docker.com/engine/tutorials/dockervolumes/


Installing from the Python Package Index (PyPI)
===============================================

The simplest way to install **iCount** is from `PyPI`_, by issuing the following command::

    pip install icount


If installing the package globally, you may need root priviliges (prefix upper command with
``sudo``).


Running within virtualenv
=========================

We recommend installing the package in `virtualenv`_. First, install the virtualenv tool and
create a virtual environment, stored in ``~/envs/icount_env``::

    pip3 install virtualenv
    virtualenv ~/envs/icount_env

Then activate the environment and install iCount into it using pip::

    source ~/envs/icount_env/bin/activate
    pip install icount

To use iCount, make sure that the proper virtualenv is loaded::

    source ~/envs/icount_env/bin/activate

Afterwards you can import iCount::

    python
    >>> import iCount

Or use its command line interface::

    iCount -h

.. _`virtualenv`:
    https://virtualenv.pypa.io/en/stable/

.. _`PyPI`:
    https://pypi.python.org/pypi


Installing from source
======================

If you wish to install form source, follow instructions in section
:doc:`Contributing <contributing>`.


