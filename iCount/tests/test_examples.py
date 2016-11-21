# pylint: disable=missing-docstring, protected-access

import os
import unittest
import tempfile
import warnings

import iCount


class TestExamplesScriptsInstall(unittest.TestCase):

    def setUp(self):
        self.tempdir = tempfile.mkdtemp()
        self.examples_dir = os.path.join(self.tempdir, 'examples')
        warnings.simplefilter("ignore", ResourceWarning)

    def test_examples(self):
        iCount.examples.run(out_dir=self.tempdir)
        # check if two scripts are present in subfolder examples
        self.assertTrue(
            os.path.isfile(os.path.join(self.examples_dir, 'hnRNPC.sh'))
        )
        self.assertTrue(
            os.path.isfile(os.path.join(self.examples_dir, 'hnRNPC_reduced.sh'))
        )

    def tearDown(self):
        files = os.listdir(self.examples_dir)
        for fn in files:
            os.remove(os.path.join(self.examples_dir, fn))

        os.rmdir(self.examples_dir)
        os.rmdir(self.tempdir)

if __name__ == '__main__':
    unittest.main()
