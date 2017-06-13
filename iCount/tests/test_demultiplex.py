# pylint: disable=missing-docstring, protected-access
import unittest
import warnings

import iCount
from iCount import demultiplex
from iCount.tests.utils import get_temp_dir, make_sequence, make_quality_scores, \
    get_temp_file_name, make_list_from_file


class TestDemultiplex(unittest.TestCase):

    def setUp(self):
        self.dir = get_temp_dir()
        warnings.simplefilter("ignore", ResourceWarning)

    def test_run_fail(self):
        message = r'Output directory does not exist. Make sure it does.'
        with self.assertRaisesRegex(FileNotFoundError, message):
            demultiplex.run('reads.fq', 'adapter', ['barcodes'], out_dir='/unexisting/dir')

    def test_run_ok_with_adapter(self):
        adapter = 'GGAACC'
        barcodes = ['NNAAAN', 'NNACTN']
        data = [
            ['@header/1', 'GGAAAG' + make_sequence(40) + adapter, '+',
             make_quality_scores(50, min_chr=33, max_chr=74) + '!J'],
            ['@header2 blah', 'TTCCTT' + make_sequence(40) + adapter, '+',
             make_quality_scores(50, min_chr=33, max_chr=74) + '!J'],
            ['@header3', 'TTGGGT' + make_sequence(40) + adapter, '+',
             make_quality_scores(50, min_chr=33, max_chr=74) + '!J'],
        ]
        fq_fname = get_temp_file_name(extension='fq')
        fq_file = iCount.files.fastq.FastqFile(fq_fname, 'wt')
        for line in data:
            fq_file.write(iCount.files.fastq.FastqEntry(*line))
        fq_file.close()

        demultiplex.run(fq_fname, adapter, barcodes, mismatches=1, out_dir=self.dir)

        fq1_list = make_list_from_file('{}/demux_{}.fastq.gz'.format(self.dir, barcodes[0]))
        expected1 = [
            ['@header:rbc:GGG/1'],
            [data[0][1][6:-6]],
            ['+'],
            [data[0][3][6:-6]],
        ]
        self.assertEqual(fq1_list, expected1)

        fq2_list = make_list_from_file('{}/demux_{}.fastq.gz'.format(self.dir, barcodes[1]))
        expected2 = [
            ['@header2:rbc:TTT'],
            [data[1][1][6:-6]],
            ['+'],
            [data[1][3][6:-6]],
        ]
        self.assertEqual(fq2_list, expected2)

        fq3_list = make_list_from_file('{}/demux_{}.fastq.gz'.format(self.dir, 'nomatch'))
        expected3 = [
            ['@header3'],
            [data[2][1]],
            ['+'],
            [data[2][3]],
        ]
        self.assertEqual(fq3_list, expected3)

    def test_run_ok_no_adapter(self):
        barcodes = ['NNAAAN', 'NNACTN']
        data = [
            ['@header/1', 'GGAAAG' + make_sequence(40), '+',
             make_quality_scores(50, min_chr=65, max_chr=73) + '!J'],
            ['@header2 blah', 'TTCCTT' + make_sequence(40), '+',
             make_quality_scores(50, min_chr=65, max_chr=73) + '!J'],
            ['@header3', 'TTGGGT' + make_sequence(40), '+',
             make_quality_scores(50, min_chr=65, max_chr=73) + '!J'],
        ]
        fq_fname = get_temp_file_name(extension='fq')
        fq_file = iCount.files.fastq.FastqFile(fq_fname, 'wt')
        for line in data:
            fq_file.write(iCount.files.fastq.FastqEntry(*line))
        fq_file.close()

        demultiplex.run(fq_fname, None, barcodes, mismatches=1, out_dir=self.dir)

        fq1_list = make_list_from_file('{}/demux_{}.fastq.gz'.format(self.dir, barcodes[0]))
        expected1 = [
            ['@header:rbc:GGG/1'],
            [data[0][1][6:]],
            ['+'],
            [data[0][3][6:]],
        ]
        self.assertEqual(fq1_list, expected1)

        fq2_list = make_list_from_file('{}/demux_{}.fastq.gz'.format(self.dir, barcodes[1]))
        expected2 = [
            ['@header2:rbc:TTT'],
            [data[1][1][6:]],
            ['+'],
            [data[1][3][6:]],
        ]
        self.assertEqual(fq2_list, expected2)

        fq3_list = make_list_from_file('{}/demux_{}.fastq.gz'.format(self.dir, 'nomatch'))
        expected3 = [
            ['@header3'],
            [data[2][1]],
            ['+'],
            [data[2][3]],
        ]
        self.assertEqual(fq3_list, expected3)


class TestExtract(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", ResourceWarning)

    def test_extract_ok_1(self):
        # No mismatch, first barcode.
        barcodes = ['NNAAAN', 'NNACTN', 'NNACGN']
        adapter = 'CCCCCC'
        data = ['@header1', 'GGAAAG' + adapter + make_sequence(40), '+',
                make_quality_scores(50) + '!J']
        fq_fname = get_temp_file_name(extension='fq')
        fq_file = iCount.files.fastq.FastqFile(fq_fname, 'wt')
        fq_file.write(iCount.files.fastq.FastqEntry(*data))
        fq_file.close()

        for read, exp_id, randomer in demultiplex._extract(fq_fname, barcodes, mismatches=1):
            self.assertEqual(exp_id, 0)
            self.assertEqual(randomer, 'GGG')
            self.assertEqual(read.id, data[0])
            self.assertEqual(read.seq, data[1][6:])
            self.assertEqual(read.plus, '+')
            self.assertEqual(read.qual, data[3][6:])

    def test_extract_ok_2(self):
        # One mismatch, second barcode.
        barcodes = ['NNAAAN', 'NNCCTN', 'NNACGN']
        adapter = 'CCCCCC'
        data = ['@header1', 'TTACTT' + adapter + make_sequence(40), '+',
                make_quality_scores(50) + '!J']
        fq_fname = get_temp_file_name(extension='fq')
        fq_file = iCount.files.fastq.FastqFile(fq_fname, 'wt')
        fq_file.write(iCount.files.fastq.FastqEntry(*data))
        fq_file.close()

        for read, exp_id, randomer in demultiplex._extract(fq_fname, barcodes, mismatches=1):
            self.assertEqual(exp_id, 1)
            self.assertEqual(randomer, 'TTT')
            self.assertEqual(read.id, data[0])
            self.assertEqual(read.seq, data[1][6:])
            self.assertEqual(read.plus, '+')
            self.assertEqual(read.qual, data[3][6:])

    def test_extract_mismatch(self):
        # To many mismatches
        barcodes = ['NNAAAN', 'NNCCTN', 'NNACGN']
        adapter = 'CCCCCC'
        data = ['@header1', 'TTACTT' + adapter + make_sequence(40), '+',
                make_quality_scores(50) + '!J']
        fq_fname = get_temp_file_name(extension='fq')
        fq_file = iCount.files.fastq.FastqFile(fq_fname, 'wt')
        fq_file.write(iCount.files.fastq.FastqEntry(*data))
        fq_file.close()

        for read, exp_id, randomer in demultiplex._extract(fq_fname, barcodes, mismatches=0):
            self.assertEqual(exp_id, -1)
            self.assertEqual(randomer, '')
            self.assertEqual(read.id, data[0])
            self.assertEqual(read.seq, data[1])
            self.assertEqual(read.plus, '+')
            self.assertEqual(read.qual, data[3])


class TestMakeP2n2i(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", ResourceWarning)

    def test_make_p2n2i(self):
        barcodes = [
            'NNAAAN',
            'NNACTN',
            'NNACGN',
        ]

        result = demultiplex._make_p2n2i(barcodes)
        expected = {
            2: {
                'A': {0, 1, 2}
            },
            3: {
                'A': {0},
                'C': {1, 2},
            },
            4: {
                'A': {0},
                'T': {1},
                'G': {2},
            },
        }
        self.assertEqual(result, expected)


if __name__ == '__main__':
    unittest.main()
