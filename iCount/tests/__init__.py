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
