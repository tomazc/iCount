"""TODO."""
# pylint: disable=missing-docstring, protected-access

import os

import iCount

FUNCTIONAL_TEST_FOLDER = os.path.join(iCount.TMP_ROOT, 'tests/functional/')

if not os.path.exists(FUNCTIONAL_TEST_FOLDER):
    print("Test folder does not exist. Will create it at: "
          "{0:s}".format(FUNCTIONAL_TEST_FOLDER))
    os.makedirs(FUNCTIONAL_TEST_FOLDER)
