************
Contributing
************


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
CLI, etc.) can be read in :doc:`CLI interface section <auto_cli>`

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
Errors shlud be raised in a descriptive manner - for example if iCount calls an
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
it in :doc:`CLI interface section <auto_cli>`.


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
