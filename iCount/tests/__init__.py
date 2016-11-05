"""
Tests for iCount python package and corresponding CLI.

There are two types of tests: *unit* and *functional* (regression) tests.
Continious development testing supported on GitHub is also explained in the
corresponding section.


Unit tests
----------

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
----------------

Are located in subdirectory iCount/tests/functional. They should be executed
manually, only at times one wishes to check that current implementation is
compatible with earlier versions of iCount. The results are typically stored for
future reference. These test may take significant amount of time to complete
and are not meant to be run on daily basis.


Continuous development testing
------------------------------

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
"""
import os
import unittest

import iCount


def suite(loader=None, pattern='test*.py'):
    test_dir = os.path.dirname(__file__)
    if loader is None:
        loader = unittest.TestLoader()
    if pattern is None:
        pattern = 'test*.py'
    top_level_dir = os.path.dirname(os.path.dirname(iCount.__file__))
    all_tests = [
        loader.discover(test_dir, pattern, top_level_dir),
    ]

    return unittest.TestSuite(all_tests)


def load_tests(loader, tests, pattern):
    return suite(loader, pattern)


if __name__ == '__main__':
    unittest.main(defaultTest='suite')
