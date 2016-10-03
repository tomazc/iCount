import os
import shutil
import unittest
import warnings

from iCount import demultiplex
from iCount.tests.utils import make_fastq_file, make_list_from_file, get_temp_dir


class TestDemultiplex(unittest.TestCase):

    def setUp(self):
        self.dir = get_temp_dir()
        self.adapter = 'CCCCCCCCC'
        self.barcodes = [
            'NNNGGTTNN',
            'NNNTTGTNN',
            'NNNCAATNN',
            'NNNACCTNN',
            'NNNGGCGNN',
        ]
        self.reads = make_fastq_file(barcodes=self.barcodes, adapter=self.adapter)
        warnings.simplefilter("ignore", ResourceWarning)

    def test_run_fail(self):
        message = r'Output directory does not exist. Make sure it does.'
        with self.assertRaisesRegex(FileNotFoundError, message):
            filenames = demultiplex.run(
                self.reads, self.adapter, 1, self.barcodes, outdir='/unexisting/dir')

    def test_run_ok(self):
        expected = ['{}/demux_{}.fastq.gz'.format(self.dir, b) for b in self.barcodes]
        expected.extend([self.dir + '/demux_nomatch.fastq.gz'])

        # Without adapter
        filenames = demultiplex.run(
            self.reads, None, 1, self.barcodes, outdir=self.dir)
        self.assertEqual(sorted(filenames), sorted(expected))

        # With adapter
        filenames = demultiplex.run(
            self.reads, self.adapter, 1, self.barcodes, outdir=self.dir)
        self.assertEqual(sorted(filenames), sorted(expected))


if __name__ == '__main__':
    unittest.main()
