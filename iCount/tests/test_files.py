import os
import gzip
import unittest
import tempfile
import warnings

import iCount


class TestFilesTemp(unittest.TestCase):

    def setUp(self):
        self.tempdir = tempfile.mkdtemp()
        warnings.simplefilter("ignore", ResourceWarning)

    def test_uncompressed(self):
        fn_in = os.path.join(self.tempdir, 'uncompressed.txt')
        # create file
        test_text = 'empty file'
        with open(fn_in, 'wt') as f:
            f.write(test_text)
        fn_out = iCount.files.decompress_to_tempfile(fn_in)
        self.assertEqual(fn_in, fn_out)
        # content should be same
        with open(fn_out, 'rt') as f:
            read_text = f.read()
        self.assertEqual(read_text, test_text)

    def test_compressed(self):
        fn_in = os.path.join(self.tempdir, 'compressed.txt.gz')
        # create file
        test_text = 'empty file'
        with gzip.open(fn_in, 'wt') as f:
            f.write(test_text)
        fn_out = iCount.files.decompress_to_tempfile(fn_in)
        self.assertNotEqual(fn_in, fn_out)
        # content should be same
        with open(fn_out, 'rt') as f:
            read_text = f.read()
        self.assertEqual(read_text, test_text)

    def tearDown(self):
        files = os.listdir(self.tempdir)
        for file_ in files:
            os.remove(os.path.join(self.tempdir, file_))
        os.rmdir(self.tempdir)

if __name__ == '__main__':
    unittest.main()
