import os
import unittest
import warnings

import pysam
from numpy import random

from iCount.mapping import xlsites

from iCount.tests.utils import get_temp_file_name, make_bam_file


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
        cur_vals = {'A': [1, 2, 3], 'B': [4, 5, 0]}
        to_add = {'A': [1, 0, 1], 'C': [42, 2, 3]}

        expected = {'A': [2, 2, 4], 'B': [4, 5, 0], 'C': [42, 2, 3]}
        xlsites._update(cur_vals, to_add)
        self.assertEqual(expected, cur_vals)


class TestMergeSimilarRandomers(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", ResourceWarning)

    def test_merge(self):
        # hit1, hit2, should be a 4-tuple, but for this test it is ok as is
        by_bc = {'AAAAA': ['hit1', 'hit2'],
                 'AAAAG': ['hit3'],
                 'TTTTN': ['hit42'],
                 'GGGGG': ['hit4', 'hit5'],
                 'NAAAN': ['hit6']}

        xlsites._merge_similar_randomers(by_bc, randomer_mismatches=2)
        expected = {'GGGGG': ['hit4', 'hit5'],
                    'AAAAA': ['hit1', 'hit2', 'hit6', 'hit3'],
                    'TTTTN': ['hit42']}
        self.assertEqual(by_bc, expected)

    @unittest.skip
    def test_todo(self):
        """
        This test fails since we currently do not support the feature of finding
        the most similar match. Rather than that, we accept the first match even
        though better one is availiable.
        """

        # hit1, hit2, should be a 4-tuple, but for this test it is ok as is
        by_bc = {'AAAAA': ['hit1'],
                 'AACGG': ['hit2', 'hit3'],
                 'NAAAN': ['hit4']}
        xlsites._merge_similar_randomers(by_bc, randomer_mismatches=2)
        expected = {'AACGG': ['hit2', 'hit3'],
                    'AAAAA': ['hit1', 'hit4']}
        self.assertEqual(by_bc, expected)


class TestCollapse(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", ResourceWarning)

    def test_start_1(self):
        xlink_pos = 1
        report_by = 'start'
        by_bc = {
            'AAAAA': [
                # (middle_pos, end_pos, read_len, num_mapped)
                (5, 10, 10, 1),
                (5, 11, 10, 1),
                (5, 12, 80, 8)],
            'CCCCC': [
                (5, 10, 10, 1),
                (5, 10, 90, 9)]}

        # Multimax = 1
        result1 = xlsites._collapse(xlink_pos, by_bc, report_by, multimax=1)
        expected1 = {1: (2.0, 3)}
        self.assertEqual(result1, expected1)

        # Multimax = 10
        result2 = xlsites._collapse(xlink_pos, by_bc, report_by, multimax=10)
        expected2 = {1: (0.5, 5)}
        self.assertEqual(result2, expected2)

    def test_middle(self):
        xlink_pos = 1
        report_by = 'middle'
        by_bc = {
            'AAAAA': [
                # (middle_pos, end_pos, read_len, num_mapped)
                (5, 10, 10, 1),
                (5, 11, 10, 1),
                (10, 12, 80, 8)],
            'CCCCC': [
                (5, 10, 10, 1),
                (20, 30, 10, 10)]}

        # Multimax = 1
        result1 = xlsites._collapse(xlink_pos, by_bc, report_by, multimax=1)
        expected1 = {5: (2.0, 3)}
        self.assertEqual(result1, expected1)

        # Multimax = 10
        result2 = xlsites._collapse(xlink_pos, by_bc, report_by, multimax=10)
        expected2 = {5: (0.7, 3), 10: (0.1, 1), 20: (0.05, 1)}
        self.assertEqual(result2, expected2)

    def test_end(self):
        xlink_pos = 1
        report_by = 'end'
        by_bc = {
            'AAAAA': [
                # (middle_pos, end_pos, read_len, num_mapped)
                (1, 5, 10, 1),
                (1, 5, 10, 1),
                (1, 10, 80, 8)],
            'CCCCC': [
                (1, 5, 10, 1),
                (1, 20, 10, 10)]}

        # Multimax = 1
        result1 = xlsites._collapse(xlink_pos, by_bc, report_by, multimax=1)
        expected1 = {5: (2.0, 3)}
        self.assertEqual(result1, expected1)

        # Multimax = 10
        result2 = xlsites._collapse(xlink_pos, by_bc, report_by, multimax=10)
        expected2 = {5: (0.7, 3), 10: (0.1, 1), 20: (0.05, 1)}
        self.assertEqual(result2, expected2)


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
                ('name3:rbc:CCCC:', 16, 0, 100, 20, [(0, 100)], {'NH': 1}),
                # Bad name - has ':' char but randomer is 'ABC' which is invalid
                # Length is divisible by 2 on positive strand
                # chr2 (RNAME=1):
                ('name4:ABC', 0, 1, 300, 20, [(0, 200)], {'NH': 11}),
                # Bad name - has ':' char but randomer is 'ACG' which is OK
                # chr2 (RNAME=1):
                ('name4:ACG', 0, 1, 300, 20, [(0, 200)], {'NH': 11}),
                # Low quality(MAPQ=3)
                ('name5', 0, 1, 300, 3, [(0, 200)], {'NH': 13})]}

        self.data_no_NH = {
            'chromosomes': [('chr1', 3000)],
            'segments': [
                # No NH tag is set
                ('name5', 0, 0, 300, 20, [(0, 200)], {})]}
        warnings.simplefilter("ignore", ResourceWarning)

    def test_run_simple(self):
        bam_fname = make_bam_file(self.data)
        unique_fname = get_temp_file_name()
        multi_fname = get_temp_file_name()

        result = xlsites.run(bam_fname, unique_fname, multi_fname, mapq_th=5)

        # All records:
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

    def test_error_open_bamfile(self):
        """
        Provide onyl file with no content - error shoud be raised.
        """
        bam_fname = get_temp_file_name()
        unique_fname = get_temp_file_name()
        multi_fname = get_temp_file_name()

        message = r"Error opening BAM file: .*"
        with self.assertRaisesRegex(ValueError, message):
            result = xlsites.run(bam_fname, unique_fname, multi_fname)

    def test_no_NH_tag(self):
        bam_fname = make_bam_file(self.data_no_NH)
        unique_fname = get_temp_file_name()
        multi_fname = get_temp_file_name()

        message = r'"NH" tag not set for record: .*'
        with self.assertRaisesRegex(ValueError, message):
            result = xlsites.run(bam_fname, unique_fname, multi_fname)


if __name__ == '__main__':
    unittest.main()
