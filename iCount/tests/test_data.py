"""
Downloads test data files from the public iCount server.
"""
# pylint: disable=missing-docstring, protected-access

import os
import urllib.request
import iCount

TEST_FOLDER = os.path.join(iCount.TMP_ROOT, 'tests')
if not os.path.exists(TEST_FOLDER):
    print("Test folder does not exist. Will create it at: "
          "{0:s}".format(TEST_FOLDER))
    os.makedirs(TEST_FOLDER)

TEST_FOLDER = os.path.join(TEST_FOLDER, 'data')
if not os.path.exists(TEST_FOLDER):
    print("Test data folder does not exist. Will create it at: "
          "{0:s}".format(TEST_FOLDER))
    os.makedirs(TEST_FOLDER)


def get_update(files_to_get):
    print('Get or update files needed for tests.')
    for fn_save, url in files_to_get:
        save_as = os.path.join(TEST_FOLDER, fn_save)

        print('downloading: {:s}'.format(url))
        if not os.path.isfile(save_as):
            urllib.request.urlretrieve(url, save_as)
            print('    done, stored as: {:s}'.format(save_as))
        else:
            print('    stored as: {:s}'.format(save_as))
    print('All files ready.')
