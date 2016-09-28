import unittest
from unittest import mock

import iCount

import iCount.externals.cutadapt as cutadapt
import iCount.externals.star as star

from iCount.tests.utils import make_fasta_file, make_fastq_file, get_temp_dir, \
    get_temp_file_name, make_file_from_list


class TestCutadapt(unittest.TestCase):

    def setUp(self):
        self.adapter = 'AAAATTTTCCCCGGGG'
        self.reads = make_fastq_file(adapter=self.adapter, num_sequences=100,
                                     out_file=get_temp_file_name(extension='fastq'))
        self.tmp = get_temp_file_name(extension='fastq')

    def test_get_version_ok(self):
        version = cutadapt.get_version()
        self.assertRegex(version, r'\d\.\d+')

    @mock.patch('iCount.externals.cutadapt.subprocess')
    def test_get_version_fail(self, subprocess_mock):
        subprocess_mock.check_output.side_effect = Exception()
        version = cutadapt.get_version()
        self.assertIsNone(version)

    def test_run(self):
        return_code = cutadapt.run(self.reads, self.tmp, self.adapter, qual_base=64,
                                   qual_trim=30, minimum_length=70)
        self.assertEqual(return_code, 0)


class TestStar(unittest.TestCase):

    def setUp(self):
        self.dir = get_temp_dir()
        self.index_dir = get_temp_dir()
        self.genome = make_fasta_file(num_sequences=2, seq_len=1000)
        self.reads = make_fastq_file(genome=self.genome)
        self.annotation = make_file_from_list([
            ['1', '.', 'gene', '10', '20', '.', '+', '.',
             'gene_id "A";'],
            ['1', '.', 'transcript', '10', '20', '.', '+', '.',
             'gene_id "A"; transcript_id "AA";'],
            ['1', '.', 'exon', '10', '20', '.', '+', '.',
             'gene_id "A"; transcript_id "AA"; exon_number "1";'],
        ])

    def test_get_version_ok(self):
        version = star.get_version()
        # Version example: STAR_2.5.0a
        regex = r'STAR_\d\.[\d\w]+'
        self.assertRegex(version, regex)

    @mock.patch('iCount.externals.star.subprocess')
    def test_get_version_fail(self, subprocess_mock):
        subprocess_mock.check_output.side_effect = Exception()
        version = star.get_version()
        self.assertIsNone(version)

    def test_build_index_bad_outdir(self):

        message = r'Output directory does not exist. Make sure it does.'
        with self.assertRaisesRegex(FileNotFoundError, message):
            star.build_index(self.genome, '/unexisting/outdir')

    def test_build_index(self):
        # No annotation
        return_code1 = star.build_index(self.genome, self.index_dir, overhang=100,
                                        overhang_min=8, threads=1, quiet=True)
        # With annotation
        return_code2 = star.build_index(self.genome, self.index_dir, annotation=self.annotation,
                                        overhang=100, overhang_min=8, threads=1, quiet=True)

        self.assertEqual(return_code1, 0)
        self.assertEqual(return_code2, 0)

    def test_map_reads_bad_genomedir(self):

        message = r'Directory with genome index does not exist. Make sure it does.'
        with self.assertRaisesRegex(FileNotFoundError, message):
            star.map_reads(self.reads, '/unexisting/genomedir', self.dir)

    def test_map_reads_bad_outdir(self):

        message = r'Output directory does not exist. Make sure it does.'
        with self.assertRaisesRegex(FileNotFoundError, message):
            star.map_reads(self.reads, self.dir, '/unexisting/outdir')

    def test_map_reads(self):
        # First: make index:
        star.build_index(self.genome, self.index_dir, quiet=True)

        # No annotation
        return_code1 = star.map_reads(self.reads, self.index_dir, self.dir, quiet=True)

        # With annotation:
        return_code2 = star.map_reads(
            self.reads, self.index_dir, self.dir, annotation=self.annotation,
            multimax=10, mismatches=2, threads=1, quiet=True)

        self.assertEqual(return_code1, 0)
        self.assertEqual(return_code2, 0)


if __name__ == '__main__':
    unittest.main()
