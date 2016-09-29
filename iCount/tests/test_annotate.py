import os
import unittest
import tempfile

import pybedtools

from iCount.analysis import annotate
from iCount.tests.utils import make_file_from_list, make_list_from_file


def template(cross_links, annotation, subtype='biotype',
             excluded_types=None):
    """
    Utility function for testing iCount.analysis.annotate

    Instead of input files, accept the file content in form of lists and create
    temporary files from them on the fly. This avoids the problem of having a
    bunch of multiple small files or one large file (which would violate the
    idea of test isolation).

    For example of how to use this function check any test that uses it.

    Parameters
    ----------
    cross_links : list
        List representation of cross-links file.
    annotation : list
        List representation of annotation file.

    Returns
    -------
    list
        List representation of output file of analysis.annotate().

    """
    cross_links_file = make_file_from_list(cross_links)
    annotation_file = make_file_from_list(annotation)
    out_file = tempfile.NamedTemporaryFile(delete=False).name
    return make_list_from_file(annotate.annotate_cross_links(
        annotation_file, cross_links_file, out_file, subtype=subtype,
        excluded_types=excluded_types), fields_separator='\t')


class TestAnnotateCrossLinks(unittest.TestCase):

    def setUp(self):
        self.annotation_file = os.path.normpath(os.path.join(
            os.path.dirname(__file__), "files/annotation.get_regions.gtf"))
        self.cross_links_file = os.path.normpath(os.path.join(
            os.path.dirname(__file__), "files/cross_links.bed"))

    def test_diff_chromosome_naming(self):
        """
        Check that exception is raised if chromosome naming is inconsistent.
        """
        cross_links = [
            ['chr1', '15', '16', '.', '5', '+'],
        ]
        annotation = [
            ['1', '.', 'CDS', '10', '20', '.', '+', '.', 'biotype "ABC";'],
        ]

        message = r"No intersections found. This may be caused by .*"
        with self.assertRaisesRegex(ValueError, message):
            template(cross_links, annotation)

    def test_basic(self):
        """Detect single cross_link and annotate it."""
        cross_links = [
            ['1', '15', '16', '.', '5', '+'],
        ]
        annotation = [
            ['1', '.', 'CDS', '10', '20', '.', '+', '.', 'biotype "ABC";'],
        ]
        expected = [
            ['1', '15', '16', 'CDS ABC', '5', '+'],
        ]
        self.assertEqual(template(cross_links, annotation), expected)

    def test_multiple_annotations(self):
        """Single cross_link intersects with multiple annotations."""
        cross_links = [
            ['1', '15', '16', '.', '5', '+'],
        ]
        annotation = [
            ['1', '.', 'CDS', '10', '20', '.', '+', '.', 'biotype "A";'],
            ['1', '.', 'CDS', '10', '20', '.', '+', '.', 'biotype "A";'],
            ['1', '.', 'CDS', '10', '20', '.', '+', '.', 'biotype "B";'],
            ['1', '.', 'CDS', '10', '20', '.', '-', '.', 'biotype "C";'],
            ['1', '.', 'CDS', '12', '18', '.', '+', '.', 'biotype "A";'],
            ['1', '.', 'CDS', '30', '40', '.', '+', '.', 'biotype "D";'],
            ['1', '.', 'ncRNA', '10', '20', '.', '+', '.', 'biotype "A";'],
        ]
        expected = [
            ['1', '15', '16', 'CDS A, CDS B, ncRNA A', '5', '+'],
        ]
        self.assertEqual(template(cross_links, annotation), expected)

    def test_shuffled_contents(self):
        """Result should be sorted, even if inputs are not."""
        cross_links = [
            ['1', '16', '17', '.', '5', '+'],
            ['1', '14', '15', '.', '5', '+'],
            ['1', '15', '16', '.', '5', '+'],
        ]
        annotation = [
            ['1', '.', 'CDS', '10', '20', '.', '+', '.', 'biotype "A";'],
        ]
        expected = [
            ['1', '14', '15', 'CDS A', '5', '+'],
            ['1', '15', '16', 'CDS A', '5', '+'],
            ['1', '16', '17', 'CDS A', '5', '+'],
        ]
        self.assertEqual(template(cross_links, annotation), expected)

    def test_strand_specific(self):
        """Two cross-links with same coordinate, but different strand."""
        cross_links = [
            ['1', '15', '16', '.', '5', '+'],
            ['1', '15', '16', '.', '5', '-'],
        ]
        annotation = [
            ['1', '.', 'intron', '10', '20', '.', '+', '.', 'biotype "A";'],
            ['1', '.', 'intron', '10', '20', '.', '+', '.', 'biotype "B";'],
            ['1', '.', 'ncRNA', '10', '20', '.', '+', '.', 'biotype "A";'],
            ['1', '.', 'ncRNA', '10', '20', '.', '-', '.', 'biotype "B";'],
        ]
        expected = [
            ['1', '15', '16', 'intron A, intron B, ncRNA A', '5', '+'],
            ['1', '15', '16', 'ncRNA B', '5', '-'],
        ]
        self.assertEqual(template(cross_links, annotation), expected)

    def test_subtype_param1(self):
        """Subtype parameter can be empty: type is just 3rd column."""
        cross_links = [
            ['1', '5', '6', '.', '1', '+'],
        ]
        annotation = [
            ['1', '.', 'intron', '1', '10', '.', '+', '.', 'biotype "A";'],
            ['1', '.', 'intron', '2', '10', '.', '+', '.', 'biotype "B";'],
            ['1', '.', 'ncRNA', '3', '10', '.', '+', '.', 'biotype "C";'],
        ]
        expected = [
            ['1', '5', '6', 'intron, ncRNA', '1', '+'],
        ]

        self.assertEqual(template(
            cross_links, annotation, subtype=None), expected)

    def test_subtype_param2(self):
        """Subtype can have any value - not just the default biotype."""
        cross_links = [
            ['1', '5', '6', '.', '1', '+'],
        ]
        annotation = [
            ['1', '.', 'intron', '1', '10', '.', '+', '.', 'attr42 "A";'],
            ['1', '.', 'intron', '2', '10', '.', '+', '.', 'attr42 "B";'],
            ['1', '.', 'ncRNA', '3', '10', '.', '+', '.', 'attr42 "C";'],
        ]
        expected = [
            ['1', '5', '6', 'intron A, intron B, ncRNA C', '1', '+'],
        ]

        self.assertEqual(template(
            cross_links, annotation, subtype='attr42'), expected)

    def test_excluded_types(self):
        """Exclude annotation intervals by 3rd column value."""
        cross_links = [
            ['1', '5', '6', '.', '1', '+'],
        ]
        annotation = [
            ['1', '.', 'intron', '1', '10', '.', '+', '.', 'biotype "A";'],
            ['1', '.', 'intron', '2', '10', '.', '+', '.', 'biotype "B";'],
            ['1', '.', 'ncRNA', '3', '10', '.', '+', '.', 'biotype "C";'],
        ]
        expected = [
            ['1', '5', '6', 'ncRNA C', '1', '+']]

        self.assertEqual(template(
            cross_links, annotation, excluded_types=['intron']), expected)


if __name__ == '__main__':
    unittest.main()
