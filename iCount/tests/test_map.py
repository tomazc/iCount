import unittest

from iCount.mapping.mapindex import run as run_mapindex
from iCount.mapping.map import run
from iCount.tests.utils import make_fasta_file, make_fastq_file, get_temp_dir


class TestMap(unittest.TestCase):

    def setUp(self):
        self.genome = make_fasta_file(num_sequences=2, seq_len=1000)
        self.reads = make_fastq_file(genome=self.genome)
        self.dir = get_temp_dir()
        self.index_dir = get_temp_dir()

    def test_map_bad_genomedir(self):

        message = r'Directory with genome index does not exist. Make sure it does.'
        with self.assertRaisesRegex(FileNotFoundError, message):
            run(self.genome, '/unexisting/genomedir', self.dir)

    def test_map_bad_outdir(self):

        message = r'Output directory does not exist. Make sure it does.'
        with self.assertRaisesRegex(FileNotFoundError, message):
            run(self.genome, self.dir, '/unexisting/outdir')

    def test_mapindex_ok(self):
        # First make index:
        run_mapindex(self.genome, self.index_dir)

        # Run the actual mapping:
        return_code = run(self.reads, self.index_dir, self.dir)
        self.assertEqual(return_code, 0)


if __name__ == '__main__':
    unittest.main()
