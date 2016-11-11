************
Contributing
************


Installation for development
============================

Fork the main `iCount git repository`_.

If you don’t have Git installed on your system, follow `these instructions`_.

Clone your fork (replace <username> with your GitHub account name) and change
directory::

    git clone https://github.com/<username>/iCount.git
    cd iCount

Prepare iCount for development::

    pip install -e .[tests]

.. note::

    We strongly recommend to install iCount in `virtualenv`_ to create an
    isolated Python environment.


.. _`iCount git repository`:
    https://github.com/tomazc/iCount

.. _`these instructions`:
    https://git-scm.com/book/en/v2/Getting-Started-Installing-Git

.. _`virtualenv`:
    https://virtualenv.pypa.io/en/stable/


Python module and CLI living together
=====================================

iCount package offers a complete set of Python function to analyse iCLIP data,
*as well* as their command line counterparts. Effectively, they (should) offer
excactly the same functionality. However, the tools are in continious
development and pipelines will change with time... To avoid differences in
functionality as well as the duplication of code, tests and documentation, the
command line interface (CLI) is generated *automatically* from Python functions.
For this to work, functions (as well as their docstrings and modules where they
are defined) should be written by the conventions described in the following
sections.

More about automatic CLI generation (how it works, how to expose function in
CLI, etc.) can be read in :obj:`Automated CLI creation <iCount.cli>`.

Conventions
===========

Syntax
------

* Line length is 99 characters plus newline character.
* Use ``'`` instead of ``"`` character, if possible.

Parameters naming
-----------------

If the variable is one of the *commonly used names*, such as cross-link file
name or mismatches, the name should be one of those. For a list of these names,
check the TODO: link.

Errors
------
Errors should be raised in a descriptive manner - for example if iCount calls an
external tool that return a non-zero exit status, an error should be
raised.

Docstrings
----------

We follow `Numpy`_ format specification for docstrings.

All docstrings of unittest test functions should have three-braces (""") in
separate lines (or else, docstrings are printed to stdout during test
execution). `This`_ explains why:

Correct way::

      def test_something(self):
          """
          The docstring.
          """

Wrong way::

    def test_something(self):
        """The docstring."""

Commit messages
---------------

If the core change that commit is introducing originates from module
``the_module``, commit message should be::

    the_module: Commit message

If commit is somehow connected with an issue, commit message should reference
issue number::

    the_module: Commit message, relates to #42, resolves #14


Logging logic
=============

Logging is used to report the progress of program execution to user. iCount can
be used as Python module or as it's corresponding CLI.

When using iCount as Python module all logging is turned OFF by defult, so no
messages will be printed about analysis execution. However, if user wishes,
logging can be easily confiugured with functions defined in
:file:`iCount.logger.py`.

If using the CLI, logging to stdout with INFO level is set by defult. This can be configured for each
command by using the appropriate command line arguents. Read more about
it in :obj:`Automated CLI creation <iCount.cli>`.


Tests
=====

.. automodule:: iCount.tests


Naming proposition - TODO
=========================

The iCount pipeline looks like so: TODO: image

The elements appearing in it are files and analysis:

    For each file (type):
        * format specification = FASTA, FASTQ, GTF, BED6, special, custom...
        * one-letter-CLI-abreviation: '-r'
        * full CLI name: '--reads'
        * naming convention: example for genome file: species.release[.chr1_chr2_chrMT].fa[.gz]
        * is output of which analysis: ?? (should be clear form image?)
        * can be input for which analysis: ?? (should be clear form image?)

    The analysis in it are:
        * analysis name (single word?)
        * execution time (O(n^2), O(reads#, genome_size)) some typical estimation
        * inputs, outputs (should be clear form image?)


.. _`Numpy`:
    http://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html

.. _`This`:
    http://stackoverflow.com/questions/12962772/how-to-stop-python-unittest-from-printing-test-docstring


Preparing a release
===================

Pull the latest master to be released::

    git checkout master
    git pull

Use the utility script ``docs/changelog.sh`` to list all commits that were made from last tagged
release::

    docs/changelog.sh > to_edit.rst

Edit ``to_edit.rst`` and incorporate a description of most important changes into ``docs/source/revisions.rst``.

Modify ``setup.py`` and set::

    ISRELEASED = True

Commit changes::

    git add docs/source/revisions.rst setup.py
    git commit -m "Release <VERSION>"

Use current value of VERSION in ``setup.py`` to tag the commit::

    git tag <VERSION>

Decide how to bump version (to some new value NEWVERSION) and modify ``setup.py``::

    VERSION=<NEWVERSION>
    ISRELEASED = False

Commit changes::

    git add setup.py
    git commit -m "Bump version to <NEWVERSION>"

Push changes to master::

    git push origin --tags


Pushing to pypi
===============

Pull a tagged version::

    git checkout <VERSION>


Create source distribution::

    python setup.py sdist

Upload to pypi using `twine`_::

    twine upload dist/iCount-<VERSION>*


.. _`twine`:
    https://pypi.python.org/pypi/twine
