************
Contributing
************

Thanks for taking the time to contribute to iCount!

Please submit contributions in accordance with the flow explained in the
`GitHub Guides`_.

.. _`GitHub Guides`:
    https://guides.github.com/


Installation for development
============================

Fork the main `iCount git repository`_.

If you don’t have Git installed on your system, follow `these instructions`_.

Clone your fork (replace <username> with your GitHub account name) and change
directory::

    git clone https://github.com/<username>/iCount.git
    cd iCount

Prepare iCount for development::

    pip install -e .[test]

.. note::

    We strongly recommend to install iCount in `virtualenv`_ to create an
    isolated Python environment.


.. _`iCount git repository`:
    https://github.com/tomazc/iCount

.. _`these instructions`:
    https://git-scm.com/book/en/v2/Getting-Started-Installing-Git

.. _`virtualenv`:
    https://virtualenv.pypa.io/en/stable/



.. _`installing-docker`:

Installing Docker
=================

When working with docker, make sure that the docker-machine has enough memory to run STAR and
associated programs, *e.g.*, at least 32 GB RAM and 45 GB disk::

    docker-machine create -d virtualbox --virtualbox-memory 32768 --virtualbox-disk-size "46080" default

.. note::
    Other options for `VirtualBox`_ are described `here`_.

.. _`VirtualBox`:
    https://www.virtualbox.org/

.. _`here`:
    https://docs.docker.com/machine/drivers/virtualbox/


Building the documentation
==========================

Issue these commands to build the documentation::

    pip install -e .[docs]
    cd docs
    make html
    open build/html/index.html

.. note::
    Calling makefile ensures that the reference manual for the CLI is updated. It will envoke
    this command when needed::

        iCount man > source/ref_CLI.txt

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

* Line length is 119 characters plus newline character.
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

Commit and PR messages
----------------------

If the core change that commit is introducing originates from module
``the_module``, commit message should be::

    the_module: Commit message

If commit is somehow connected with an issue, commit message should reference
issue number::

    the_module: Commit message, relates to #42, resolves #14

References to PR greatly help when preparing a description of changes in the change log ``docs/source/revisions.rst``.


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

Tests for iCount python package and corresponding CLI.
------------------------------------------------------

There are two types of tests: *unit* and *functional* (regression) tests.
Continuous development testing supported on GitHub is also explained in the
corresponding section.


Unit tests
^^^^^^^^^^

Unit tests are located in top testing directory (iCount/tests/). They follow the
standardy philosophy of unit tests (isolation), although this is sometimes
violated when testing functions connect to external resources. This is only the
case, when the main task of the function under test is to retrieve a resource
from the web. Still, all tests should pass in no more than a couple of minutes.

Test can be run by standard unittest call::

    # Run all tests:
    python -m unittest iCount
    # Run only a specific test:
    python -m unittest iCount.test_file_name.Test_class.test_name

Alternatively, all test files can be called also like python script::

    python test_file.py


Functional tests
^^^^^^^^^^^^^^^^

Are located in subdirectory iCount/tests/functional. They should be executed
manually, only at times one wishes to check that current implementation is
compatible with earlier versions of iCount. The results are typically stored for
future reference. These test may take significant amount of time to complete
and are not meant to be run on daily basis.


Continuous development testing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

iCount is project in continious development and therefore has central repository
on GitHub. There, automatic testing is performed by:

    * Travis - unit, code and documentation stype testing
    * Coverage - enforcing sufficient covergae
    * Codacy - enforcing code quality
    * Scrutinizer - enforcing documentation quality

Tests on Travis are executed by tox. To avoid making multiple pull requests and
waiting for Travis to do the checks, tox can run on local machine::

    # Run all tox tests
    tox
    # Run only a single environment (for example unittests in Python 3.5)
    tox -e py35

Tox enables to test code against different versions of Python and also perform
code-style checks. This is achieved by different *environments*. Currently, we
test three of them:

    * py35 (Python 3.5)
    * py34 (Python 3.4)
    * linters (wraps pylint, PEP8/pycodestyle and PEP257/pydocstyle)

For more info check ``tox.ini`` and ``.travis.yml`` in packyge root directory.


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

First check that ``twine`` is installed::

    pip install twine

    Requirement already satisfied: twine in ...

Pull the latest master to be released::

    git checkout master
    git pull

Use the utility script ``docs/changelog.sh`` to list all commits that were made from last tagged
release::

    docs/changelog.sh > to_edit.rst

Edit ``to_edit.rst`` and incorporate a description of most important changes into ``docs/source/revisions.rst``.
Use syntax from `releases`_ package.

.. _`releases`:
    http://releases.readthedocs.io/en/latest/concepts.html#issue-and-release-types

Remove ``.dev`` from project's version in ``iCount/__about__.py`` file.

Clean ``build`` directory::

    python setup.py clean -a

Remove previous distributions in ``dist`` directory::

    rm dist/*

Remove previous ``egg-info`` directory::

    rm -r *.egg-info

Commit changes::

    git add docs/source/revisions.rst iCount/__about__.py
    git commit -m "Release <version>"

Test the new version with Tox_::

    tox -r

Create source distribution::

    python setup.py sdist

Build wheel::

    python setup.py bdist_wheel

Upload distribution to PyPI_::

    twine upload dist/*

.. note::
    It is advisable to test upload onto the test server https://test.pypi.org/legacy/ first::

        twine upload --repository-url https://test.pypi.org/legacy/ dist/*

    Afterwards, test the installation in a clean python environment::

        docker build -t icounttestinstall -f Dockerfile_test .
        docker run -t -i icounttestinstall

        pip install -i https://test.pypi.org/pypi \
        --extra-index-url https://pypi.python.org/pypi iCount

Tag the new version::

    git tag <version>

Push changes to main repository::

    git push origin master <version>

Decide how to bump version (to some new value <new-version>) and modify 
``iCount/__about__.py``. Don't forget to add suffix ``.dev``::

    __version__=<new-version>.dev

.. note::

    Use `Semantic versioning`_.

Commit changes::

    git add iCount/__about__.py
    git commit -m "Bump version to <new-version>"


.. _`twine`:
    https://pypi.python.org/pypi/twine

.. _Semantic versioning:
    https://packaging.python.org/en/latest/distributing/#semantic-versioning-preferred

.. _Tox:
    http://tox.testrun.org/
     
.. _PyPi:
    https://pypi.python.org/

