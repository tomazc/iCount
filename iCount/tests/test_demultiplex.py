# pylint: disable=missing-docstring, protected-access

import os
import unittest
import warnings

import iCount
from iCount import demultiplex
from iCount.files.fastq import FastqEntry, FastqFile
from iCount.tests.utils import get_temp_dir, make_sequence, make_quality_scores, \
    get_temp_file_name, make_list_from_file


class TestCreateIndex(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", ResourceWarning)

    def test_create_index_simple(self):
        barcodes5 = [
            'NNAAAN',
            'NNACTN',
            'NNACGN',
        ]
        barcodes = demultiplex.prepare_barcodes(barcodes5, None)

        index = demultiplex.create_index(barcodes)
        expected = {
            2: {
                'A': {barcodes5[0], barcodes5[1], barcodes5[2]}
            },
            3: {
                'A': {barcodes5[0]},
                'C': {barcodes5[1], barcodes5[2]},
            },
            4: {
                'A': {barcodes5[0]},
                'T': {barcodes5[1]},
                'G': {barcodes5[2]},
            },
        }
        self.assertEqual(index, expected)

    def test_barcode_size_diff(self):
        barcodes5 = [
            'NATAN',
            'NNNCAGN',
        ]
        barcodes = demultiplex.prepare_barcodes(barcodes5, None)

        index = demultiplex.create_index(barcodes)
        expected = {
            1: {
                'A': {barcodes5[0]}
            },
            2: {
                'T': {barcodes5[0]}
            },
            3: {
                'A': {barcodes5[0]},
                'C': {barcodes5[1]},
            },
            4: {
                'A': {barcodes5[1]},
            },
            5: {
                'G': {barcodes5[1]},
            },
        }
        self.assertEqual(index, expected)

    def test_barcode3(self):
        barcodes5 = [
            'NATN',
            'NNCAN',
            'NNCAN',
        ]
        barcodes3 = [
            'NNTAG',
            'NNGT',
            '.',
        ]

        barcodes = demultiplex.prepare_barcodes(barcodes5, barcodes3)
        index = demultiplex.create_index(barcodes)
        self.assertEqual(index, {
            1: {
                'A': {barcodes5[0]}
            },
            2: {
                'T': {barcodes5[0]},
                'C': {barcodes5[1]},
            },
            3: {
                'A': {barcodes5[1]},
            },
        })

        index3 = demultiplex.create_index(barcodes['NATN']['barcodes3'])
        self.assertEqual(index3, {
            -3: {
                'T': {barcodes3[0]},
            },
            -2: {
                'A': {barcodes3[0]},
            },
            -1: {
                'G': {barcodes3[0]},
            },
        })

        index3 = demultiplex.create_index(barcodes['NNCAN']['barcodes3'])
        self.assertEqual(index3, {
            -2: {
                'G': {barcodes3[1]},
            },
            -1: {
                'T': {barcodes3[1]},
            },
        })


class TestExtract(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", ResourceWarning)

        self.barcodes5 = [
            'NNAAAN',
            'NNACTN',
            'NNACGANN',
        ]
        self.barcodes3 = [
            '.',
            '.',
            'NNGGG',
            'NNAAA',
        ]

    def create_fq_file(self, entries):
        fname = get_temp_file_name(extension='fq')

        fq_file = FastqFile(fname, 'wt')
        for entry in entries:
            fq_file.write(entry)
        fq_file.close()
        return fname

    def test_extract_ok_1(self):
        # No mismatch, first barcode.
        barcodes = demultiplex.prepare_barcodes(self.barcodes5, None)
        entry = FastqEntry(
            '@header1',
            'GGAAAG' + make_sequence(40),
            '+',
            make_quality_scores(46),
        )
        fq_fname = self.create_fq_file([entry])

        for read, wbrc, randomer in demultiplex._extract(fq_fname, barcodes, mismatches=1, minimum_length=15):
            self.assertEqual(wbrc, self.barcodes5[0])
            self.assertEqual(randomer, 'GGG')
            self.assertEqual(read.id, entry.id)
            self.assertEqual(read.seq, entry.seq[6:])
            self.assertEqual(read.plus, entry.plus)
            self.assertEqual(read.qual, entry.qual[6:])

    def test_extract_ok_2(self):
        # One mismatch, second barcode.
        barcodes = demultiplex.prepare_barcodes(self.barcodes5, None)
        entry = FastqEntry(
            '@header1',
            'TTAGTT' + make_sequence(40),
            '+',
            make_quality_scores(46),
        )
        fq_fname = self.create_fq_file([entry])

        for read, wbrc, randomer in demultiplex._extract(fq_fname, barcodes, mismatches=1, minimum_length=15):
            self.assertEqual(wbrc, self.barcodes5[1])
            self.assertEqual(randomer, 'TTT')
            self.assertEqual(read.id, entry.id)
            self.assertEqual(read.seq, entry.seq[6:])
            self.assertEqual(read.plus, entry.plus)
            self.assertEqual(read.qual, entry.qual[6:])

    def test_extract_mismatch(self):
        # To many mismatches
        barcodes = demultiplex.prepare_barcodes(self.barcodes5, None)
        entry = FastqEntry(
            '@header1',
            'TTTTTT' + make_sequence(40),
            '+',
            make_quality_scores(46),
        )
        fq_fname = self.create_fq_file([entry])

        for read, wbrc, randomer in demultiplex._extract(fq_fname, barcodes, mismatches=0, minimum_length=15):
            self.assertEqual(wbrc, 'nomatch')
            self.assertEqual(randomer, '')
            self.assertEqual(read.id, entry.id)
            self.assertEqual(read.seq, entry.seq)
            self.assertEqual(read.plus, entry.plus)
            self.assertEqual(read.qual, entry.qual)

    def test_barcode_size_diff(self):
        # One mismatch, second barcode.
        barcodes = demultiplex.prepare_barcodes(self.barcodes5, None)
        entry1 = FastqEntry(
            '@header1',
            'TTAGAT' + make_sequence(40),
            '+',
            make_quality_scores(44),
        )
        entry2 = FastqEntry(
            '@header2',
            'GGACGAGG' + make_sequence(40),
            '+',
            make_quality_scores(48),
        )
        fq_fname = self.create_fq_file([entry1, entry2])

        handle = demultiplex._extract(fq_fname, barcodes, mismatches=1, minimum_length=15)
        read1, wbrc, randomer1 = next(handle)
        self.assertEqual(wbrc, self.barcodes5[0])
        self.assertEqual(randomer1, 'TTT')
        self.assertEqual(read1.id, entry1.id)
        self.assertEqual(read1.seq, entry1.seq[6:])
        self.assertEqual(read1.plus, entry1.plus)
        self.assertEqual(read1.qual, entry1.qual[6:])
        read2, wbrc, randomer2 = next(handle)
        self.assertEqual(wbrc, self.barcodes5[2])
        self.assertEqual(randomer2, 'GGGG')
        self.assertEqual(read2.id, entry2.id)
        self.assertEqual(read2.seq, entry2.seq[8:])
        self.assertEqual(read2.plus, entry2.plus)
        self.assertEqual(read2.qual, entry2.qual[8:])

    def test_barcode3(self):
        # One mismatch, second barcode.
        barcodes = demultiplex.prepare_barcodes(self.barcodes5 + [self.barcodes5[2]], self.barcodes3)
        barcodes = barcodes[self.barcodes5[2]]['barcodes3']
        entry1 = FastqEntry(
            '@header1',
            make_sequence(40) + 'TTGGG',
            '+',
            make_quality_scores(45),
        )
        entry2 = FastqEntry(
            '@header2',
            make_sequence(40) + 'TCAAA',
            '+',
            make_quality_scores(45),
        )
        fq_fname = self.create_fq_file([entry1, entry2])

        handle = demultiplex._extract(fq_fname, barcodes, mismatches=1, minimum_length=15)
        read1, wbrc, randomer1 = next(handle)
        self.assertEqual(wbrc, self.barcodes3[2])
        self.assertEqual(randomer1, 'TT')
        self.assertEqual(read1.id, entry1.id)
        self.assertEqual(read1.seq, entry1.seq[:-5])
        self.assertEqual(read1.plus, entry1.plus)
        self.assertEqual(read1.qual, entry1.qual[:-5])
        read2, wbrc, randomer2 = next(handle)
        self.assertEqual(wbrc, self.barcodes3[3])
        self.assertEqual(randomer2, 'TC')
        self.assertEqual(read2.id, entry2.id)
        self.assertEqual(read2.seq, entry2.seq[:-5])
        self.assertEqual(read2.plus, entry2.plus)
        self.assertEqual(read2.qual, entry2.qual[:-5])


class TestRun(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", ResourceWarning)

        self.dir = get_temp_dir()
        self.adapter = 'AAAAAAAAAA'
        self.barcodes5 = [
            'NNAAAN',
            'NGGGN',
            'NGGGN',
        ]
        self.barcodes3 = [
            '.',
            'NNGGG',
            'NCCC',
        ]
        # Header: early version Illumina header
        # Barcodes: exact match to the barcode set #1
        self.entry1 = FastqEntry(
            '@header1/1',
            'GGAAAG' + make_sequence(40) + self.adapter,
            '+',
            make_quality_scores(56),
        )
        # Header: contains id and description
        # Barcodes: one mismatch on 5' end for barcode set #2
        self.entry2 = FastqEntry(
            '@header2 blah',
            'AGGTA' + make_sequence(40) + 'AAGGG' + self.adapter,
            '+',
            make_quality_scores(60),
        )
        # Header: simple header
        # Barcodes: one mismatch on 3' end for barcode set #3
        self.entry3 = FastqEntry(
            '@header3',
            'TGGGT' + make_sequence(40) + 'TACC' + self.adapter,
            '+',
            make_quality_scores(59),
        )

        self.fq_fname = get_temp_file_name(extension='fq')
        self.fq_file = iCount.files.fastq.FastqFile(self.fq_fname, 'wt')
        for entry in [self.entry1, self.entry2, self.entry3]:
            self.fq_file.write(entry)
        self.fq_file.close()

    def test_no_out_dir(self):
        message = r'Output directory does not exist. Make sure it does.'
        with self.assertRaisesRegex(FileNotFoundError, message):
            demultiplex.run('reads.fq', 'adapter', ['barcodes'], out_dir='/unexisting/dir')

    def test_only_barcode5_0_mismatch(self):
        # Only barcode5, one mismatch
        demultiplex.run(self.fq_fname, self.adapter, self.barcodes5[:2], mismatches=0, out_dir=self.dir)

        demux_file = 'demux_{}.fastq.gz'.format(self.barcodes5[0])
        fq_list = make_list_from_file(os.path.join(self.dir, demux_file))
        self.assertEqual(fq_list[0], ['@header1:rbc:GGG/1'])
        self.assertEqual(fq_list[1], [self.entry1.seq[6:-10]])
        self.assertEqual(fq_list[3], [self.entry1.qual[6:-10]])

        demux_file = 'demux_{}.fastq.gz'.format(self.barcodes5[1])
        fq_list = make_list_from_file(os.path.join(self.dir, demux_file))
        self.assertEqual(fq_list[0], ['@header3:rbc:TT'])
        self.assertEqual(fq_list[1], [self.entry3.seq[5:-10]])
        self.assertEqual(fq_list[3], [self.entry3.qual[5:-10]])

        demux_file = 'demux_{}.fastq.gz'.format('nomatch5')
        fq_list = make_list_from_file(os.path.join(self.dir, demux_file))
        self.assertEqual(fq_list[0], ['@header2'])
        self.assertEqual(fq_list[1], [self.entry2.seq])
        self.assertEqual(fq_list[3], [self.entry2.qual])

    def test_only_barcode5_1_mismatch(self):
        # Only barcode5, one mismatch
        demultiplex.run(self.fq_fname, self.adapter, self.barcodes5[:2], mismatches=1, out_dir=self.dir)

        demux_file = 'demux_{}.fastq.gz'.format(self.barcodes5[0])
        fq_list = make_list_from_file(os.path.join(self.dir, demux_file))
        self.assertEqual(fq_list[0], ['@header1:rbc:GGG/1'])
        self.assertEqual(fq_list[1], [self.entry1.seq[6:-10]])
        self.assertEqual(fq_list[3], [self.entry1.qual[6:-10]])

        demux_file = 'demux_{}.fastq.gz'.format(self.barcodes5[1])
        fq_list = make_list_from_file(os.path.join(self.dir, demux_file))
        self.assertEqual(fq_list[0], ['@header2:rbc:AA'])
        self.assertEqual(fq_list[1], [self.entry2.seq[5:-10]])
        self.assertEqual(fq_list[3], [self.entry2.qual[5:-10]])
        self.assertEqual(fq_list[4], ['@header3:rbc:TT'])
        self.assertEqual(fq_list[5], [self.entry3.seq[5:-10]])
        self.assertEqual(fq_list[7], [self.entry3.qual[5:-10]])

        demux_file = 'demux_{}.fastq.gz'.format('nomatch5')
        fq_list = make_list_from_file(os.path.join(self.dir, demux_file))
        self.assertEqual(fq_list, [])

    def test_both_barcodes_1_mismatch(self):
        # Both: barcode5 & barcode3, one mismatch
        demultiplex.run(self.fq_fname, self.adapter, self.barcodes5,
                        barcodes3=self.barcodes3, mismatches=1, out_dir=self.dir)

        demux_file = 'demux_{}.fastq.gz'.format(self.barcodes5[0])
        fq_list = make_list_from_file(os.path.join(self.dir, demux_file))
        self.assertEqual(fq_list[0], ['@header1:rbc:GGG/1'])
        self.assertEqual(fq_list[1], [self.entry1.seq[6:-10]])
        self.assertEqual(fq_list[3], [self.entry1.qual[6:-10]])

        demux_file = 'demux_{}_{}.fastq.gz'.format(self.barcodes5[1], self.barcodes3[1])
        fq_list = make_list_from_file(os.path.join(self.dir, demux_file))
        self.assertEqual(len(fq_list), 4)
        self.assertEqual(fq_list[0], ['@header2:rbc:AAAA'])
        self.assertEqual(fq_list[1], [self.entry2.seq[5:-15]])
        self.assertEqual(fq_list[3], [self.entry2.qual[5:-15]])

        demux_file = 'demux_{}_{}.fastq.gz'.format(self.barcodes5[2], self.barcodes3[2])
        fq_list = make_list_from_file(os.path.join(self.dir, demux_file))
        self.assertEqual(fq_list, [])

        demux_file = 'demux_{}.fastq.gz'.format('nomatch5')
        fq_list = make_list_from_file(os.path.join(self.dir, demux_file))
        self.assertEqual(fq_list, [])

        demux_file = 'demux_{}_{}.fastq.gz'.format(self.barcodes5[1], 'nomatch')
        fq_list = make_list_from_file(os.path.join(self.dir, demux_file))
        self.assertEqual(len(fq_list), 4)
        self.assertEqual(fq_list[0], ['@header3:rbc:TT'])
        self.assertEqual(fq_list[1], [self.entry3.seq[5:-10]])
        self.assertEqual(fq_list[3], [self.entry3.qual[5:-10]])

        demux_file = 'no_adapter_found_{}.fastq.gz'.format(self.barcodes5[1])
        fq_list = make_list_from_file(os.path.join(self.dir, demux_file))
        self.assertEqual(fq_list, [])


if __name__ == '__main__':
    unittest.main()
