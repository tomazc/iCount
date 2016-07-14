"""
Downloads test data files from the public iCount server.


"""

import os
import urllib.request
import iCount

test_folder = os.path.join(iCount.storage_root, 'tests')
if not os.path.exists(test_folder):
    print("Test folder does not exist. Will create it at: "
          "{0:s}".format(test_folder))
    os.makedirs(test_folder)

test_data_folder = os.path.join(test_folder, 'data')
if not os.path.exists(test_data_folder):
    print("Test data folder does not exist. Will create it at: "
          "{0:s}".format(test_data_folder))
    os.makedirs(test_data_folder)


def get_update(files_to_get):
    print('Get or update files needed for tests.')
    for fn_save, url in files_to_get:
        save_as = os.path.join(test_data_folder, fn_save)

        print('downloading: {:s}'.format(url))
        if not os.path.isfile(save_as):
            urllib.request.urlretrieve(url, save_as)
            print('    done, stored as: {:s}'.format(save_as))
        else:
            print('    stored as: {:s}'.format(save_as))
    print('All files ready.')
