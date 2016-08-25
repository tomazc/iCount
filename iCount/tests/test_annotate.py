import os
import unittest
import tempfile

import pybedtools

from iCount.analysis import annotate


def make_file_from_list(data):
    """Return path to file with the content from `data` (list of lists)"""
    tfile = tempfile.NamedTemporaryFile(delete=False)
    pybedtools.BedTool(pybedtools.create_interval_from_list(list_)
                       for list_ in data).saveas(tfile.name)
    tfile.close()
    return os.path.abspath(tfile.name)


def make_list_from_file(fname, fields_separator=None):
    """Read file to list of lists"""
    data = []
    with open(fname) as file_:
        for line in file_:
            data.append(line.strip().split(fields_separator))
    return data


def test_template(cross_links, annotation):
    """
    Utility function for testing iCount.analysis.annotate

    Instead of input files, accept the file content in form of lists and create
    temporary files from them on the fly. This avoids the problem of having a
    bunch of multiple small files or one large file (which would violate the
    idea of test isolation).

    For example of how to use this function check any test that uses it.

    :param list cross_links: list representation of cross-links file
    :param list annotation: list representation of annotation file
    :return: list representation of output file of analysis.annotate()
    :rtype: list
    """
    cross_links_file = make_file_from_list(cross_links)
    annotation_file = make_file_from_list(annotation)
    out_file = tempfile.NamedTemporaryFile(delete=False).name
    return make_list_from_file(annotate.annotate_cross_links(
        annotation_file, cross_links_file, out_file), fields_separator='\t')


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
            ['chr1', '15', '16', '.', '5', '+']]
        annotation = [
            ['1', '.', 'CDS', '10', '20', '.', '+', '.', 'biotype "ABC";']]

        message = r"No intersections found. This may be caused by .*"
        with self.assertRaisesRegex(ValueError, message):
            test_template(cross_links, annotation)

    def test_basic(self):
        """Detect single cross_link and annotate it."""
        cross_links = [
            ['1', '15', '16', '.', '5', '+']]
        annotation = [
            ['1', '.', 'CDS', '10', '20', '.', '+', '.', 'biotype "ABC";']]
        expected = [
            ['1', '15', '16', 'CDS ABC', '5', '+']]
        self.assertEqual(test_template(cross_links, annotation), expected)

    def test_multiple_annotations(self):
        """Single cross_link intersects with multiple annotations."""
        cross_links = [
            ['1', '15', '16', '.', '5', '+']]
        annotation = [
            ['1', '.', 'CDS', '10', '20', '.', '+', '.', 'biotype "A";'],
            ['1', '.', 'CDS', '10', '20', '.', '+', '.', 'biotype "A";'],
            ['1', '.', 'CDS', '10', '20', '.', '+', '.', 'biotype "B";'],
            ['1', '.', 'CDS', '10', '20', '.', '-', '.', 'biotype "C";'],
            ['1', '.', 'CDS', '12', '18', '.', '+', '.', 'biotype "A";'],
            ['1', '.', 'CDS', '30', '40', '.', '+', '.', 'biotype "D";'],
            ['1', '.', 'ncRNA', '10', '20', '.', '+', '.', 'biotype "A";']]
        expected = [
            ['1', '15', '16', 'CDS A, CDS B, ncRNA A', '5', '+']]
        self.assertEqual(test_template(cross_links, annotation), expected)

    def test_shuffled_contents(self):
        """Result should be sorted, even if inputs are not."""
        cross_links = [
            ['1', '16', '17', '.', '5', '+'],
            ['1', '14', '15', '.', '5', '+'],
            ['1', '15', '16', '.', '5', '+']]
        annotation = [
            ['1', '.', 'CDS', '10', '20', '.', '+', '.', 'biotype "A";']]
        expected = [
            ['1', '14', '15', 'CDS A', '5', '+'],
            ['1', '15', '16', 'CDS A', '5', '+'],
            ['1', '16', '17', 'CDS A', '5', '+']]
        self.assertEqual(test_template(cross_links, annotation), expected)

    def test_strand_specific(self):
        """Two cross-links with same coordinate, but different strand."""
        cross_links = [
            ['1', '15', '16', '.', '5', '+'],
            ['1', '15', '16', '.', '5', '-']]
        annotation = [
            ['1', '.', 'intron', '10', '20', '.', '+', '.', 'biotype "A";'],
            ['1', '.', 'intron', '10', '20', '.', '+', '.', 'biotype "B";'],
            ['1', '.', 'ncRNA', '10', '20', '.', '+', '.', 'biotype "A";'],
            ['1', '.', 'ncRNA', '10', '20', '.', '-', '.', 'biotype "B";']]
        expected = [
            ['1', '15', '16', 'intron A, intron B, ncRNA A', '5', '+'],
            ['1', '15', '16', 'ncRNA B', '5', '-']]
        self.assertEqual(test_template(cross_links, annotation), expected)


if __name__ == '__main__':
    unittest.main()
