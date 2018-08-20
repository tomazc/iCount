# pylint: disable=missing-docstring, protected-access

import os
import gzip
import unittest
import tempfile
import warnings

import iCount
from iCount.tests.utils import get_temp_file_name, make_file_from_list, make_list_from_file


class TestFilesTemp(unittest.TestCase):

    def setUp(self):
        self.tempdir = tempfile.mkdtemp()
        warnings.simplefilter("ignore", ResourceWarning)

    def test_uncompressed(self):
        fn_in = os.path.join(self.tempdir, 'uncompressed.txt')
        # create file
        test_text = 'empty file'
        with open(fn_in, 'wt') as file_:
            file_.write(test_text)
        fn_out = iCount.files.decompress_to_tempfile(fn_in)
        self.assertEqual(fn_in, fn_out)
        # content should be same
        with open(fn_out, 'rt') as file_:
            read_text = file_.read()
        self.assertEqual(read_text, test_text)

    def test_compressed(self):
        fn_in = os.path.join(self.tempdir, 'compressed.txt.gz')
        # create file
        test_text = 'empty file'
        with gzip.open(fn_in, 'wt') as file_:
            file_.write(test_text)
        fn_out = iCount.files.decompress_to_tempfile(fn_in)
        self.assertNotEqual(fn_in, fn_out)
        # content should be same
        with open(fn_out, 'rt') as file_:
            read_text = file_.read()
        self.assertEqual(read_text, test_text)

    def tearDown(self):
        files = os.listdir(self.tempdir)
        for file_ in files:
            os.remove(os.path.join(self.tempdir, file_))
        os.rmdir(self.tempdir)


class TestFilesFastq(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", ResourceWarning)

    def test_get_qual_encoding_s(self):
        data = [['@header'], ['AAAAA'], ['+'], ['!FFFI']]
        fq_file = make_file_from_list(data, bedtool=False, extension='.fq')
        self.assertEqual('S', iCount.files.fastq.get_qual_encoding(fq_file))

    def test_get_qual_encoding_x(self):
        data = [['@header'], ['AAAAA'], ['+'], [';FFFh']]
        fq_file = make_file_from_list(data, bedtool=False, extension='.fq')
        self.assertEqual('X', iCount.files.fastq.get_qual_encoding(fq_file))

    def test_get_qual_encoding_i(self):
        data = [['@header'], ['AAAAA'], ['+'], ['@FFFh']]
        fq_file = make_file_from_list(data, bedtool=False, extension='.fq')
        self.assertEqual('I', iCount.files.fastq.get_qual_encoding(fq_file))

    def test_get_qual_encoding_j(self):
        data = [['@header'], ['AAAAA'], ['+'], ['BFFFh']]
        fq_file = make_file_from_list(data, bedtool=False, extension='.fq')
        self.assertEqual('J', iCount.files.fastq.get_qual_encoding(fq_file))

    def test_get_qual_encoding_l(self):
        data = [['@header'], ['AAAAA'], ['+'], ['!FFFJ']]
        fq_file = make_file_from_list(data, bedtool=False, extension='.fq')
        self.assertEqual('L', iCount.files.fastq.get_qual_encoding(fq_file))

    def test_get_qual_encoding_none(self):
        data = [['@header'], ['AAAAA'], ['+'], ['???']]
        fq_file = make_file_from_list(data, bedtool=False, extension='.fq')
        self.assertEqual(None, iCount.files.fastq.get_qual_encoding(fq_file))

    def test_get_qual_encoding_empty(self):
        data = []
        fq_file = make_file_from_list(data, bedtool=False, extension='.fq')
        self.assertEqual(None, iCount.files.fastq.get_qual_encoding(fq_file))

    def test_fastq_entry(self):
        fq_entry = iCount.files.fastq.FastqEntry('header12345', 'ACTGCTGCAT', '+', 'ABCDEFFBBA')
        self.assertEqual(fq_entry.id, 'header12345')
        self.assertEqual(fq_entry.seq, 'ACTGCTGCAT')
        self.assertEqual(fq_entry.plus, '+')
        self.assertEqual(fq_entry.qual, 'ABCDEFFBBA')

    def test_fastq_entry_1(self):
        # Illumina before 1.8
        fq_entry = iCount.files.fastq.FastqEntry(
            '@HWUSI-EAS100R:6:73:941:1973#0/1', 'AAA', '+', 'FFF')
        self.assertEqual(fq_entry.id, '@HWUSI-EAS100R:6:73:941:1973#0/1')

    def test_fastq_entry_2(self):
        # Illumina after 1.8
        fq_entry = iCount.files.fastq.FastqEntry(
            '@EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG', 'AAA', '+', 'FFF')
        self.assertEqual(fq_entry.id, '@EAS139:136:FC706VJ:2:2104:15343:197393')

    def test_fastq_entry_3(self):
        # NCBI
        fq_entry = iCount.files.fastq.FastqEntry(
            '@SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=36', 'AAA', '+', 'FFF')
        self.assertEqual(fq_entry.id, '@SRR001666.1')

    def test_fastq_file_not_exist(self):
        with self.assertRaises(FileNotFoundError):
            iCount.files.fastq.FastqFile('/nonexisting/file', 'rt')

    def test_fastq_file_read_l(self):
        data = [
            ['@EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG'],
            ['AAA'],
            ['+'],
            ['ABh'],
            ['@EAS139:136:FC706VJ:2:2104:15343:197394 1:Y:19:CTCACG'],
            ['AAAA'],
            ['+'],
            [';ABC'],
        ]
        fq_file = make_file_from_list(data, bedtool=False, extension='.fq')
        reader = iCount.files.fastq.FastqFile(fq_file, 'rt').read()
        entry1 = next(reader)
        entry2 = next(reader)
        self.assertEqual(data[0][0][:39], entry1.id)
        self.assertEqual(data[1][0], entry1.seq)
        self.assertEqual(data[2][0], entry1.plus)
        self.assertEqual(data[3][0], entry1.qual)
        self.assertEqual(data[4][0][:39], entry2.id)
        self.assertEqual(data[5][0], entry2.seq)

        with self.assertRaises(StopIteration):
            next(reader)

    def test_fastq_file_read_x(self):
        data = [
            ['@EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG'],
            ['AAA'],
            ['+'],
            ['!FI'],
            ['@EAS139:136:FC706VJ:2:2104:15343:197394 1:Y:19:CTCACG'],
            ['AAAA'],
            ['+'],
            ['!FIF'],
        ]
        fq_file = make_file_from_list(data, bedtool=False, extension='.fq')
        reader = iCount.files.fastq.FastqFile(fq_file, 'rt').read()
        entry1 = next(reader)
        entry2 = next(reader)
        self.assertEqual(data[0][0][:39], entry1.id)
        self.assertEqual(data[1][0], entry1.seq)
        self.assertEqual(data[2][0], entry1.plus)
        self.assertEqual(data[3][0], entry1.qual)
        self.assertEqual(data[4][0][:39], entry2.id)
        self.assertEqual(data[5][0], entry2.seq)

        with self.assertRaises(StopIteration):
            next(reader)

    def test_fastq_file_write(self):
        data = [
            ['@header1', 'AAA', '+', 'FFF'],
            ['@header2', 'AAAA', '+', 'FFFF'],
        ]
        fq_file_name = get_temp_file_name(extension='fq.gz')
        fq_file = iCount.files.fastq.FastqFile(fq_file_name, 'wt')
        for line in data:
            fq_file.write(iCount.files.fastq.FastqEntry(*line))
        fq_file.close()
        result = make_list_from_file(fq_file_name)
        expected = [['@header1'], ['AAA'], ['+'], ['FFF'], ['@header2'], ['AAAA'], ['+'], ['FFFF']]
        self.assertEqual(result, expected)


class TestBedGraph(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", ResourceWarning)

        bed_data = [
            ['1', '4', '5', '.', '5', '+'],
            ['1', '5', '6', '.', '1', '+'],
            ['1', '5', '6', '.', '1', '-'],
            ['2', '5', '6', '.', '3', '+'],
        ]
        self.bed = make_file_from_list(bed_data, extension='bed')
        self.bedgraph = get_temp_file_name(extension='bedgraph')

    def test_bed2bedgraph(self):
        iCount.files.bedgraph.bed2bedgraph(self.bed, self.bedgraph)
        expected = [
            ['track type=bedGraph name="User Track" description="User Supplied Track"'],
            ['1', '4', '5', '+5'],
            ['1', '5', '6', '+1'],
            ['1', '5', '6', '-1'],
            ['2', '5', '6', '+3'],
        ]
        result = make_list_from_file(self.bedgraph, fields_separator='\t')
        self.assertEqual(result, expected)

    def test_bed2bedgraph_params(self):
        """
        Test with custom parameters.
        """
        iCount.files.bedgraph.bed2bedgraph(
            self.bed,
            self.bedgraph,
            name='Sample name',
            description='A long and detailed description.',
            visibility='full',
            priority=20,
            color='256,0,0',
            alt_color='0,256,0',
            max_height_pixels='100:50:0',
        )
        expected = [
            ['track type=bedGraph name="Sample name" description="A long and detailed description."'
             ' visibility=full priority=20 color=256,0,0 altColor=0,256,0 maxHeightPixels=100:50:0'],
            ['1', '4', '5', '+5'],
            ['1', '5', '6', '+1'],
            ['1', '5', '6', '-1'],
            ['2', '5', '6', '+3'],
        ]
        result = make_list_from_file(self.bedgraph, fields_separator='\t')
        self.assertEqual(result, expected)


if __name__ == '__main__':
    unittest.main()
