import unittest

from iCount.mapping.mapindex import run
from iCount.tests.utils import make_fasta_file, get_temp_dir


class TestMapindex(unittest.TestCase):

    def setUp(self):
        self.genome = make_fasta_file(num_sequences=2, seq_len=1000)
        self.dir = get_temp_dir()

    def test_mapindex_bad_outdir(self):

        message = r'Output directory does not exist. Make sure it does.'
        with self.assertRaisesRegex(FileNotFoundError, message):
            run(self.genome, '/unexisting/outdir')

    def test_mapindex_ok(self):

        return_code = run(self.genome, self.dir)
        self.assertEqual(return_code, 0)


if __name__ == '__main__':
    unittest.main()
