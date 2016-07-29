import os

import iCount

functional_test_folder = os.path.join(iCount.tmp_root, 'tests/functional/')
if not os.path.exists(functional_test_folder):
    print("Test folder does not exist. Will create it at: "
          "{0:s}".format(functional_test_folder))
    os.makedirs(functional_test_folder)
