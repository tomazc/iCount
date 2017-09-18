# pylint: disable=missing-docstring, protected-access

import warnings
import unittest
from unittest import mock

import pybedtools

from iCount.mapping import xlsites
from iCount.tests.utils import get_temp_file_name, make_bam_file


class TestGetRandomBarcode(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", ResourceWarning)

    def test_good_name(self):
        name = xlsites._get_random_barcode('_____:rbc:AAA:____', mock.MagicMock())
        self.assertEqual(name, 'AAA')

    def test_no_rbc_key_valid_nucs(self):
        metrics = mock.MagicMock()
        name = xlsites._get_random_barcode('_____:AAA', metrics)
        self.assertEqual(name, 'AAA')

    def test_no_rbc_key_invalid_nucs(self):
        metrics = mock.MagicMock()
        metrics.invalidrandomer_recs = 0
        name = xlsites._get_random_barcode('_____:AAB', metrics)
        self.assertEqual(name, '')
        self.assertEqual(metrics.invalidrandomer_recs, 1)

    def test_bad_name(self):
        metrics = mock.MagicMock()
        metrics.norandomer_recs = 0
        name = xlsites._get_random_barcode('blah', metrics)
        self.assertEqual(name, '')
        self.assertEqual(metrics.norandomer_recs, 1)


class TestMatch(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", ResourceWarning)

    def test_match_basic(self):
        self.assertFalse(xlsites._match('ACGT', 'ACGG', 0))
        self.assertTrue(xlsites._match('ACGT', 'ACGG', 1))

    def test_match_unequal_len(self):
        self.assertFalse(xlsites._match('AAAAA', 'AGNN', 1))
        self.assertTrue(xlsites._match('AGNN', 'AAAAA', 2))

    def test_capital_letters(self):
        self.assertTrue(xlsites._match('AAAA', 'aaaa', 0))
        self.assertTrue(xlsites._match('AAAA', 'anna', 0))
        self.assertFalse(xlsites._match('AAAG', 'aaaa', 0))
        self.assertTrue(xlsites._match('AAAG', 'aaaa', 1))
        self.assertTrue(xlsites._match('AANG', 'aaaa', 1))

        self.assertFalse(xlsites._match('AACGG', 'NAAAN', 1))


class TestUpdate(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", ResourceWarning)

    def test_match_basic(self):
        cur_vals = {
            'A': [1, 2, 3],
            'B': [4, 5, 0],
        }
        to_add = {
            'A': [1, 0, 1],
            'C': [42, 2, 3],
        }

        expected = {
            'A': [2, 2, 4],
            'B': [4, 5, 0],
            'C': [42, 2, 3],
        }
        xlsites._update(cur_vals, to_add)
        self.assertEqual(expected, cur_vals)


class TestMergeSimilarRandomers(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", ResourceWarning)

    def test_merge(self):
        # hit1, hit2, should be a 4-tuple, but for this test it is ok as is
        by_bc = {
            'AAAAA': ['hit1', 'hit2'],
            'AAAAG': ['hit3'],
            'TTTTN': ['hit42'],
            'GGGGG': ['hit4', 'hit5'],
            'NAAAN': ['hit6'],
        }
        expected = {
            'GGGGG': ['hit4', 'hit5'],
            'AAAAA': ['hit1', 'hit2', 'hit6', 'hit3'],
            'TTTTN': ['hit42'],
        }
        xlsites._merge_similar_randomers(by_bc, mismatches=2)
        self.assertEqual(by_bc, expected)

    @unittest.skip
    def test_todo(self):
        """
        This test fails since we currently do not support the feature of finding
        the most similar match. Rather than that, we accept the first match even
        though better one is availiable.
        """

        # hit1, hit2, should be a 4-tuple, but for this test it is ok as is
        by_bc = {
            'AAAAA': ['hit1'],
            'AACGG': ['hit2', 'hit3'],
            'NAAAN': ['hit4'],
        }
        xlsites._merge_similar_randomers(by_bc, mismatches=2)
        expected = {
            'AACGG': ['hit2', 'hit3'],
            'AAAAA': ['hit1', 'hit4'],
        }
        self.assertEqual(by_bc, expected)


class TestCollapse(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", ResourceWarning)

    def test_start_1(self):
        """
        Two barcodes, one with two second_start groups.
        No multimax situation yet.
        """
        by_bc = {
            'AAAAA': [
                # (middle_pos, end_pos, read_len, num_mapped, cigar, second_start)
                (5, 10, 10, 1, 0),
                (5, 10, 10, 1, 0),
                (5, 30, 20, 1, 20),
            ],
            'CCCCC': [
                (5, 10, 10, 1, 0),
                (5, 10, 10, 1, 0),
            ],
        }
        result1 = xlsites._collapse(1, by_bc, 'start', multimax=1)
        expected1 = {1: (3.0, 5)}
        self.assertEqual(result1, expected1)

    def test_start_2(self):
        """
        Two barcodes, multimax case.
        """
        by_bc = {
            'AAAAA': [
                # (middle_pos, end_pos, read_len, num_mapped, cigar, second_start)
                (5, 10, 10, 1, 0),
                (5, 10, 10, 1, 0),
                (15, 30, 20, 2, 0),
            ],
            'CCCCC': [
                (5, 10, 10, 1, 0),
                (5, 10, 10, 100, 0),  # Excluded becouse of multimax
            ],
        }
        result1 = xlsites._collapse(1, by_bc, 'start', multimax=10)
        expected1 = {1: (1.75, 4)}
        self.assertEqual(result1, expected1)

    def test_middle(self):
        xlink_pos = 1
        report_by = 'middle'
        by_bc = {
            'AAAAA': [
                # (middle_pos, end_pos, read_len, num_mapped, cigar, second_start)
                (5, 10, 10, 1, 0),
                (5, 10, 10, 1, 0),
                (40, 80, 80, 8, 0)],
            'CCCCC': [
                (5, 10, 10, 1, 0),
                (25, 30, 10, 10, 0)]}

        # Multimax = 1
        result1 = xlsites._collapse(xlink_pos, by_bc, report_by, multimax=1)
        expected1 = {5: (2.0, 3)}
        self.assertEqual(result1, expected1)

        # Multimax = 10
        result2 = xlsites._collapse(xlink_pos, by_bc, report_by, multimax=10)
        expected2 = {5: (0.7, 3), 40: (0.1, 1), 25: (0.05, 1)}
        self.assertEqual(result2, expected2)

    def test_end(self):
        xlink_pos = 1
        report_by = 'end'
        by_bc = {
            'AAAAA': [
                # (middle_pos, end_pos, read_len, num_mapped, cigar, second_start)
                (1, 10, 10, 1, 0),
                (1, 10, 10, 1, 0),
                (1, 80, 80, 8, 0),
            ],
            'CCCCC': [
                (1, 10, 10, 1, 0),
                (15, 20, 10, 10, 0),
            ],
        }

        # Multimax = 1
        result1 = xlsites._collapse(xlink_pos, by_bc, report_by, multimax=1)
        expected1 = {10: (2.0, 3)}
        self.assertEqual(result1, expected1)

        # Multimax = 10
        result2 = xlsites._collapse(xlink_pos, by_bc, report_by, multimax=10)
        expected2 = {10: (0.7, 3), 80: (0.1, 1), 20: (0.05, 1)}
        self.assertEqual(result2, expected2)


class TestIntersectsWithAnnotaton(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", ResourceWarning)

    def test_pos_strand(self):
        segmentation = {
            'gene_id_001': {
                'tr_id_0001': [
                    mock.MagicMock(start=100),
                ],
                'gene_segment': [],
            }
        }
        self.assertTrue(
            xlsites._intersects_with_annotaton(100, segmentation, 1, '+'))
        self.assertFalse(
            xlsites._intersects_with_annotaton(101, segmentation, 1, '+'))

    def test_neg_strand(self):
        segmentation = {
            'gene_id_002': {
                'tr_id_0003': [
                    mock.MagicMock(stop=100),
                ]
            },
        }
        self.assertTrue(
            xlsites._intersects_with_annotaton(100, segmentation, 2, '-'))
        self.assertFalse(
            xlsites._intersects_with_annotaton(101, segmentation, 2, '-'))


class TestSecondStart(unittest.TestCase):

    def setUp(self):
        self.metrics = mock.MagicMock()
        self.tmp = get_temp_file_name(extension='bam.gz')
        warnings.simplefilter("ignore", ResourceWarning)

    def test_second_start_segmentation(self):
        segmentation = {
            'G001': {
                'gene_segment': [],
                'T0001': [
                    pybedtools.create_interval_from_list(
                        ['1', '.', 'exon', '100', '200', '.', '+', '.',
                         'gene_id: "G001"', 'transcript_id: "T0001"']),
                ],
            },
            'G002': {
                'gene_segment': [],
                'T0002': [
                    pybedtools.create_interval_from_list(
                        ['1', '.', 'exon', '50', '100', '.', '-', '.',
                         'gene_id: "G001"', 'transcript_id: "T0001"']),
                ],
            },
        }

        second_start, _ = xlsites._second_start(
            read=0, poss=(1, 2, 99, 100), strand='+', chrom=1,
            segmentation=segmentation, holesize_th=4)
        self.assertEqual(second_start, 99)

        second_start, _ = xlsites._second_start(
            read=0, poss=(99, 100, 199, 200), strand='-', chrom=1,
            segmentation=segmentation, holesize_th=4)
        self.assertEqual(second_start, 100)

        second_start, _ = xlsites._second_start(
            read=0, poss=(1, 2, 4, 5), strand='-', chrom=1,
            segmentation=segmentation, holesize_th=4)
        self.assertEqual(second_start, 2)

    def test_second_start_no_seg(self):
        # If hole size is lower than holesize_th, strange should be empty:
        _, is_strange = xlsites._second_start(
            read='the_read', poss=(1, 2, 5, 6), strand='+', chrom=1,
            segmentation=None, holesize_th=1)
        self.assertTrue(is_strange)

        # If hole size is lower than holesize_th, strange should be empty:
        _, is_strange = xlsites._second_start(
            read='the_read', poss=(1, 2, 5, 6), strand='+', chrom=1,
            segmentation=None, holesize_th=2)
        self.assertFalse(is_strange)


class TestProcessBamFile(unittest.TestCase):

    def setUp(self):
        self.metrics = mock.MagicMock()
        self.tmp = get_temp_file_name(extension='bam.gz')
        warnings.simplefilter("ignore", ResourceWarning)

    def test_unmapped(self):
        """
        Unmapped read (FLAG=4):
        """
        bam_fname = make_bam_file({
            'chromosomes': [('chr1', 3000)],
            'segments': [('name1', 4, 0, 0, 0, [(0, 0)], {})],
        }, rnd_seed=0)
        self.metrics.all_recs = 0
        self.metrics.notmapped_recs = 0
        self.metrics.used_recs = 0
        list(xlsites._processs_bam_file(bam_fname, self.metrics, 0, self.tmp))
        self.assertEqual(self.metrics.notmapped_recs, 1)
        self.assertEqual(self.metrics.all_recs, 1)
        self.assertEqual(self.metrics.used_recs, 0)

    def test_low_quality(self):
        """
        Unmapped read (FLAG=4):
        """
        bam_fname = make_bam_file({
            'chromosomes': [('chr1', 3000)],
            'segments': [('name1', 0, 0, 0, 3, [(0, 100)], {})],
        }, rnd_seed=0)
        self.metrics.lowmapq_recs = 0
        self.metrics.used_recs = 0
        list(xlsites._processs_bam_file(bam_fname, self.metrics, 10, self.tmp))
        self.assertEqual(self.metrics.lowmapq_recs, 1)
        self.assertEqual(self.metrics.used_recs, 0)

    def test_no_nh_tag(self):
        data_no_nh = {
            'chromosomes': [('chr1', 3000)],
            'segments': [
                # No NH tag is set
                ('name5', 0, 0, 0, 50, [(0, 100)], {})]}
        bam_fname = make_bam_file(data_no_nh, rnd_seed=0)

        message = r'"NH" tag not set for record: .*'
        with self.assertRaisesRegex(ValueError, message):
            list(xlsites._processs_bam_file(bam_fname, self.metrics, 10, self.tmp))

    def test_pos_neg_strand(self):
        """
        positive, negative strnad
        len %2 == 0  na negativnem stnadu
        and not %2
        """
        bam_fname = make_bam_file({
            'chromosomes': [('chr1', 3000)],
            'segments': [
                # (qname, flag, refname, pos, mapq, cigar, tags)
                ('_:rbc:AAA', 16, 0, 50, 255, [(0, 100)], {'NH': 1}),
                ('_:rbc:CCC', 0, 0, 50, 255, [(0, 101)], {'NH': 1}),
            ],
        }, rnd_seed=0)
        grouped = list(xlsites._processs_bam_file(bam_fname, self.metrics, 10, self.tmp))

        expected = [
            (('chr1', '+'), 0.0167, {49: {'CCC': [(100, 150, 101, 1, 0)]}}),
            (('chr1', '-'), 1.0, {150: {'AAA': [(99, 50, 100, 1, 0)]}}),
        ]
        self.assertEqual(grouped, expected)

    def test_diff_barcodes(self):
        """
        Different barcodes on same position.
        """
        bam_fname = make_bam_file({
            'chromosomes': [('chr1', 3000)],
            'segments': [
                # (qname, flag, refname, pos, mapq, cigar, tags)
                ('_:rbc:AAA', 0, 0, 50, 255, [(0, 101)], {'NH': 1}),
                ('_:rbc:CCC', 0, 0, 50, 255, [(0, 101)], {'NH': 1}),
                ('_:rbc:CCC', 0, 0, 50, 255, [(0, 101)], {'NH': 1}),
                ('_:rbc:GGG', 0, 0, 50, 255, [(0, 101)], {'NH': 1}),
            ],
        }, rnd_seed=0)
        grouped = list(xlsites._processs_bam_file(bam_fname, self.metrics, 10, self.tmp))

        expected = [
            (('chr1', '+'), 0.0167, {
                49: {
                    'AAA': [(100, 150, 101, 1, 0)],
                    'CCC': [(100, 150, 101, 1, 0), (100, 150, 101, 1, 0)],
                    'GGG': [(100, 150, 101, 1, 0)],
                }
            }),
        ]
        self.assertEqual(grouped, expected)


class TestRun(unittest.TestCase):

    def setUp(self):
        self.data = {
            'chromosomes': [('chr1', 3000), ('chr2', 2000)],
            'segments': [
                # Unmapped read (FLAG=4):
                ('name1', 4, 0, 100, 20, [(0, 100)], {'NH': 1}),
                # Name not containing ':' or ':rbc:'
                # Mappped read and + strand (FLAG=0)
                # chr1 (RNAME=0):
                ('name2', 0, 0, 100, 20, [(0, 201)], {'NH': 7}),
                # Correct name (contains ':rbc:')
                # Mappped read on - strand (FLAG=16)
                # chr1 (RNAME=0):
                ('name3:rbc:CCCC:', 16, 0, 100, 20, [(0, 50), (3, 20), (0, 50)], {'NH': 1}),
                # Bad name - has ':' char but randomer is 'ABC' which is invalid
                # Length is divisible by 2 on positive strand
                # chr2 (RNAME=1):
                ('name4:ABC', 0, 1, 300, 20, [(0, 200)], {'NH': 11}),
                # Bad name - has ':' char but randomer is 'ACG' which is OK
                # chr2 (RNAME=1):
                ('name4:ACG', 0, 1, 300, 20, [(0, 200)], {'NH': 11}),
                # Low quality(MAPQ=3)
                ('name5', 0, 1, 300, 3, [(0, 200)], {'NH': 13})]}

        warnings.simplefilter("ignore", ResourceWarning)

    def test_run_simple(self):
        bam_fname = make_bam_file(self.data, rnd_seed=0)
        unique_fname = get_temp_file_name(extension='bed.gz')
        multi_fname = get_temp_file_name(extension='bed.gz')
        strange_fname = get_temp_file_name(extension='bam')

        result = xlsites.run(
            bam_fname, unique_fname, multi_fname, strange_fname, mapq_th=5, report_progress=True)

        # pylint: disable=no-member
        self.assertEqual(result.all_recs, 6)
        # Unmapped records:
        self.assertEqual(result.notmapped_recs, 1)
        # Mapped records:
        self.assertEqual(result.mapped_recs, 5)
        # Records with too low quality:
        self.assertEqual(result.lowmapq_recs, 1)
        # Records used in analysis
        self.assertEqual(result.used_recs, 4)
        # Records with invalid randomers
        self.assertEqual(result.invalidrandomer_recs, 1)
        # Records with no randomers:
        self.assertEqual(result.norandomer_recs, 1)
        # Barcode counter:
        self.assertEqual(result.bc_cn, {'': 2, 'ACG': 1, 'CCCC': 1})
        # Strange counter:
        self.assertEqual(result.strange_recs, 1)


if __name__ == '__main__':
    unittest.main()
