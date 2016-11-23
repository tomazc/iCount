# pylint: disable=missing-docstring, protected-access

import unittest
import tempfile
import warnings

from iCount.analysis import summary
from iCount.tests.utils import make_file_from_list, make_list_from_file, get_temp_file_name


def _make_types_length(annotation, subtype='biotype', excluded_types=None):
    """
    Run function `make_types_length_file` with data from `annotation`.
    """
    annotation_file = make_file_from_list(annotation)
    out_file = get_temp_file_name()
    fai = make_file_from_list(bedtool=False, data=[
        ['1', '100'],
        ['6', '100'],
        ['20', '100'],
    ])
    result, _ = summary.make_types_length_file(
        annotation_file, fai, out_file, subtype=subtype, excluded_types=excluded_types)
    return make_list_from_file(result, fields_separator='\t')


def _make_summary_report(annotation, cross_links, chrom_lengths,
                         subtype='biotype', excluded_types=None):
    """
    Run function `make_summary_report` with input/output data as lists.
    """
    annotation_file = make_file_from_list(annotation)
    cross_links_file = make_file_from_list(cross_links)
    chrom_length_file = make_file_from_list(chrom_lengths, bedtool=False)
    out_file = tempfile.NamedTemporaryFile(delete=False).name

    return make_list_from_file(summary.make_summary_report(
        annotation_file, cross_links_file, out_file, chrom_length_file,
        subtype=subtype, excluded_types=excluded_types), fields_separator='\t')


class TestMakeTypesLengthFile(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", ResourceWarning)

    def test_merge_same_types(self):
        """
        Test that merging same type is done as expected.

        Confirm that:

            * same interval is not counted twice
            * overlapping/touching intervals are merged and only then counted
        """
        annotation = [
            ['1', '.', 'CDS', '10', '20', '.', '+', '.', 'biotype "A";'],
            ['1', '.', 'CDS', '10', '20', '.', '+', '.', 'biotype "A";'],
            ['1', '.', 'CDS', '15', '25', '.', '+', '.', 'biotype "A";'],
            ['1', '.', 'CDS', '20', '29', '.', '+', '.', 'biotype "A";']]
        expected = [
            ['CDS A', '20'],
            ['unannotated', '580'],
        ]

        self.assertEqual(expected, _make_types_length(annotation))

    def test_merge_respect_strand(self):
        """
        Test that merging is sensitive to strand.

        Intervals differing only in strand are counted separately
        """
        annotation = [
            ['1', '.', 'CDS', '10', '19', '.', '+', '.', 'biotype "A";'],
            ['1', '.', 'CDS', '10', '19', '.', '-', '.', 'biotype "A";']]
        expected = [
            ['CDS A', '20'],
            ['unannotated', '580'],
        ]

        self.assertEqual(expected, _make_types_length(annotation))

    def test_mixed_types(self):
        """
        Defect different types in same position correctly.
        """
        annotation = [
            ['1', '.', 'intron', '10', '19', '.', '+', '.', 'biotype "A";'],
            ['1', '.', 'intron', '20', '29', '.', '+', '.', 'biotype "A";'],
            ['1', '.', 'intron', '10', '19', '.', '+', '.', 'biotype "B";'],
            ['1', '.', 'ncRNA', '10', '19', '.', '+', '.', 'biotype "A";'],
            ['1', '.', 'ncRNA', '10', '19', '.', '+', '.', 'biotype "C";']]
        expected = [
            ['intron A', '20'],
            ['intron B', '10'],
            ['ncRNA A', '10'],
            ['ncRNA C', '10'],
            ['unannotated', '580'],
        ]

        self.assertEqual(expected, _make_types_length(annotation))

    def test_shuffled_input(self):
        """
        Unsorted input does not make difference.
        """
        annotation = [
            ['20', '.', 'intron', '20', '29', '.', '+', '.', 'biotype "A";'],
            ['1', '.', 'intron', '20', '29', '.', '+', '.', 'biotype "B";'],
            ['1', '.', 'intron', '10', '19', '.', '+', '.', 'biotype "A";'],
            ['6', '.', 'ncRNA', '10', '19', '.', '+', '.', 'biotype "C";']]
        expected = [
            ['intron A', '20'],
            ['intron B', '10'],
            ['ncRNA C', '10'],
            ['unannotated', '560'],
        ]

        self.assertEqual(expected, _make_types_length(annotation))

    def test_subtype_param1(self):
        """
        Subtype parameter can be empty: type is just 3rd column.
        """
        annotation = [
            ['20', '.', 'intron', '20', '29', '.', '+', '.', 'biotype "A";'],
            ['6', '.', 'ncRNA', '10', '19', '.', '+', '.', 'biotype "C";']]
        expected = [
            ['intron', '10'],
            ['ncRNA', '10'],
            ['unannotated', '580']
        ]

        self.assertEqual(expected, _make_types_length(
            annotation, subtype=None))

    def test_subtype_param2(self):
        """
        Subtype can have any value - not just the default biotype.
        """
        annotation = [
            ['20', '.', 'intron', '20', '29', '.', '+', '.', 'attr42 "A";'],
            ['6', '.', 'ncRNA', '10', '19', '.', '+', '.', 'attr42 "C";']]
        expected = [
            ['intron A', '10'],
            ['ncRNA C', '10'],
            ['unannotated', '580'],
        ]

        self.assertEqual(expected, _make_types_length(
            annotation, subtype='attr42'))

    def test_excluded_types(self):
        """
        Exclude some annotation intervals by 3rd column value.
        """
        annotation = [
            ['20', '.', 'intron', '20', '29', '.', '+', '.', 'biotype "A";'],
            ['6', '.', 'ncRNA', '10', '19', '.', '+', '.', 'biotype "C";']]
        expected = [
            ['ncRNA C', '10'],
            ['unannotated', '590'],
        ]

        self.assertEqual(expected, _make_types_length(
            annotation, excluded_types=['intron']))


class TestMakeSummaryReport(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", ResourceWarning)
        self.chrom_lengths = [['1', '100']]
        self.header = ['type', 'length', 'length %',
                       'sites #', 'sites %', 'sites enrichment',
                       'events #', 'events %', 'events enrichment']

    def test_diff_chromosome_naming(self):
        """
        Exception is raised if chromosome naming is inconsistent.
        """
        cross_links = [
            ['chr1', '15', '16', '.', '5', '+']]
        annotation = [
            ['1', '.', 'CDS', '10', '20', '.', '+', '.', 'biotype "A";']]
        message = r"No intersections found. This may be caused by .*"
        with self.assertRaisesRegex(ValueError, message):
            _make_summary_report(annotation, cross_links, self.chrom_lengths)

    def test_diff_only_strand1(self):
        """
        Same coords but diff strand and same type.
        """
        cross_links = [
            ['1', '5', '6', '.', '3', '+'],
            ['1', '5', '6', '.', '2', '-']]
        annotation = [
            ['1', '.', 'CDS', '1', '10', '.', '+', '.', 'biotype "A";'],
            ['1', '.', 'CDS', '1', '10', '.', '-', '.', 'biotype "A";']]
        expected = [
            self.header,
            ['CDS A', '20', '0.1', '2', '1.0', '10.0', '5', '1.0', '10.0'],
            ['unannotated', '180', '0.9', '0', '0.0', '0.0', '0', '0.0', '0.0'],
        ]

        out = _make_summary_report(annotation, cross_links, self.chrom_lengths)
        self.assertEqual(out, expected)

    def test_diff_only_strand2(self):
        """
        Same coords but diff strand and diff type.
        """
        cross_links = [
            ['1', '5', '6', '.', '3', '+'],
            ['1', '5', '6', '.', '1', '-']]
        annotation = [
            ['1', '.', 'CDS', '1', '10', '.', '+', '.', 'biotype "A";'],
            ['1', '.', 'CDS', '1', '10', '.', '-', '.', 'biotype "B";']]
        expected = [
            self.header,
            ['CDS A', '10', '0.05', '1', '0.5', '10.0', '3', '0.75', '15.0'],
            ['CDS B', '10', '0.05', '1', '0.5', '10.0', '1', '0.25', '5.0'],
            ['unannotated', '180', '0.9', '0', '0.0', '0.0', '0', '0.0', '0.0'],
        ]

        out = _make_summary_report(annotation, cross_links, self.chrom_lengths)
        self.assertEqual(out, expected)

    def test_many_regions(self):
        """
        Multiple annotation regions intersecting one crosslink.
        """
        cross_links = [
            ['1', '5', '6', '.', '1', '+']]
        annotation = [
            ['1', '.', 'CDS', '1', '10', '.', '+', '.', 'biotype "A";'],
            ['1', '.', 'CDS', '1', '20', '.', '+', '.', 'biotype "B";'],
            ['1', '.', 'ncRNA', '1', '10', '.', '+', '.', 'biotype "A";'],
            ['1', '.', 'ncRNA', '1', '20', '.', '+', '.', 'biotype "C";']]
        expected = [
            self.header,
            ['CDS A', '10', '0.05', '1', '0.25', '5.0', '1', '0.25', '5.0'],
            ['CDS B', '20', '0.1', '1', '0.25', '2.5', '1', '0.25', '2.5'],
            ['ncRNA A', '10', '0.05', '1', '0.25', '5.0', '1', '0.25', '5.0'],
            ['ncRNA C', '20', '0.1', '1', '0.25', '2.5', '1', '0.25', '2.5'],
            ['unannotated', '180', '0.9', '0', '0.0', '0.0', '0', '0.0', '0.0'],
        ]

        out = _make_summary_report(annotation, cross_links, self.chrom_lengths)
        self.assertEqual(out, expected)

    def test_unsorted_input(self):
        """
        Unsoreted input should make no difference.
        """
        cross_links = [
            ['1', '7', '8', '.', '2', '+'],
            ['1', '5', '6', '.', '1', '-']]
        annotation = [
            ['1', '.', 'CDS', '20', '29', '.', '+', '.', 'biotype "A";'],
            ['1', '.', 'CDS', '6', '10', '.', '+', '.', 'biotype "A";'],
            ['1', '.', 'CDS', '1', '10', '.', '-', '.', 'biotype "A";']]
        expected = [
            self.header,
            ['CDS A', '25', '0.125', '2', '1.0', '8.0', '3', '1.0', '8.0'],
            ['unannotated', '175', '0.875', '0', '0.0', '0.0', '0', '0.0', '0.0'],
        ]

        out = _make_summary_report(annotation, cross_links, self.chrom_lengths)
        self.assertEqual(out, expected)

    def test_subtype_param1(self):
        """
        Subtype parameter can be empty: type is just 3rd column.
        """
        cross_links = [
            ['1', '5', '6', '.', '1', '+']]
        annotation = [
            ['1', '.', 'CDS', '1', '10', '.', '+', '.', 'biotype "A";'],
            ['1', '.', 'CDS', '2', '10', '.', '+', '.', 'biotype "B";']]
        expected = [
            self.header,
            ['CDS', '10', '0.05', '1', '1.0', '20.0', '1', '1.0', '20.0'],
            ['unannotated', '190', '0.95', '0', '0.0', '0.0', '0', '0.0', '0.0'],
        ]

        out = _make_summary_report(annotation, cross_links,
                                   self.chrom_lengths, subtype=None)
        self.assertEqual(out, expected)

    def test_subtype_param2(self):
        """
        Subtype can have any value - not just the default biotype.
        """
        cross_links = [
            ['1', '5', '6', '.', '1', '+']]
        annotation = [
            ['1', '.', 'CDS', '1', '10', '.', '+', '.', 'attr42 "A";'],
            ['1', '.', 'CDS', '1', '10', '.', '+', '.', 'attr42 "B";'],
            ['1', '.', 'intron', '7', '10', '.', '+', '.', 'attr42 "B";']]
        expected = [
            self.header,
            ['CDS A', '10', '0.05', '1', '0.5', '10.0', '1', '0.5', '10.0'],
            ['CDS B', '10', '0.05', '1', '0.5', '10.0', '1', '0.5', '10.0'],
            ['intron B', '4', '0.02', '0', '0.0', '0.0', '0', '0.0', '0.0'],
            ['unannotated', '190', '0.95', '0', '0.0', '0.0', '0', '0.0', '0.0'],
        ]

        out = _make_summary_report(annotation, cross_links,
                                   self.chrom_lengths, subtype='attr42')
        self.assertEqual(out, expected)

    def test_excluded_types(self):
        """
        Exclude some annotation intervals by 3rd column value.
        """
        cross_links = [
            ['1', '5', '6', '.', '1', '+']]
        annotation = [
            ['1', '.', 'intron', '1', '20', '.', '+', '.', 'biotype "A";'],
            ['1', '.', 'ncRNA', '1', '20', '.', '+', '.', 'biotype "C";']]
        expected = [
            self.header,
            ['ncRNA C', '20', '0.1', '1', '1.0', '10.0', '1', '1.0', '10.0'],
            ['unannotated', '180', '0.9', '0', '0.0', '0.0', '0', '0.0', '0.0'],
        ]

        out = _make_summary_report(annotation, cross_links,
                                   self.chrom_lengths, excluded_types=['intron'])
        self.assertEqual(out, expected)


if __name__ == '__main__':
    unittest.main()
