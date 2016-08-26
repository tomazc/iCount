import os
import unittest
import tempfile

import pybedtools

from iCount.files.bed import merge_bed


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


def test_merge_bed(data):
    """
    TODO
    """
    files = []
    for file_ in data:
        files.append(make_file_from_list(file_))
    out_file = tempfile.NamedTemporaryFile(delete=False).name
    return make_list_from_file(
        merge_bed(files, out_file), fields_separator='\t')


class TestMergeBed(unittest.TestCase):

    def setUp(self):
        pass

    def test_no_inputs(self):
        """Raise if no files are given."""
        message = "At least one element expected in files list, but none found."
        with self.assertRaisesRegex(ValueError, message):
            merge_bed([], 'outfile_name')

    def test_bad_file_name(self):
        """Raise error if bad filename given as input"""
        message = r"File .* not found."
        with self.assertRaisesRegex(ValueError, message):
            merge_bed(['/bad/path/to/unexisting/file'], 'outfile_name')

    def test_merge_single_file(self):
        """Test merging of single file"""
        bed1 = [
            ['1', '4', '5', '.', '1', '+'],
            ['1', '5', '6', '.', '1', '+'],
            ['1', '5', '6', '.', '1', '-'],
            ['2', '5', '6', '.', '1', '+'],
            ['1', '5', '6', '.', '2', '+'],
            ['1', '8', '9', '.', '1', '+']]
        expected = [
            ['1', '4', '5', '.', '1', '+'],
            ['1', '5', '6', '.', '3', '+'],
            ['1', '5', '6', '.', '1', '-'],
            ['1', '8', '9', '.', '1', '+'],
            ['2', '5', '6', '.', '1', '+']]
        out = test_merge_bed([bed1])

        self.assertEqual(out, expected)

    def test_merge_more_files(self):
        """Test that merging in sone well"""
        bed1 = [
            ['1', '4', '5', '.', '1', '+'],
            ['1', '5', '6', '.', '1', '+'],
            ['1', '5', '6', '.', '1', '-']]
        bed2 = [
            ['1', '5', '6', '.', '2', '+'],
            ['2', '5', '6', '.', '1', '+'],
            ['1', '8', '9', '.', '1', '+']]
        expected = [
            ['1', '4', '5', '.', '1', '+'],
            ['1', '5', '6', '.', '3', '+'],
            ['1', '5', '6', '.', '1', '-'],
            ['1', '8', '9', '.', '1', '+'],
            ['2', '5', '6', '.', '1', '+']]
        out = test_merge_bed([bed1, bed2])

        self.assertEqual(out, expected)


if __name__ == '__main__':
    unittest.main()
