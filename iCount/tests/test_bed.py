import unittest
import tempfile
import warnings

from iCount.files.bed import merge_bed
from iCount.tests.utils import make_file_from_list, make_list_from_file


def merge_bed_wrapper(data):
    """
    TODO
    """
    files = []
    for file_ in data:
        files.append(make_file_from_list(file_))
    out_file = tempfile.NamedTemporaryFile(delete=False).name
    return make_list_from_file(
        merge_bed(out_file, files), fields_separator='\t')


class TestMergeBed(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", ResourceWarning)

    def test_no_inputs(self):
        """
        Raise if no files are given.
        """
        message = "At least one element expected in files list, but none found."
        with self.assertRaisesRegex(ValueError, message):
            merge_bed('outfile_name', [])

    def test_bad_file_name(self):
        """
        Raise error if bad filename given as input.
        """
        message = r"File .* not found."
        with self.assertRaisesRegex(ValueError, message):
            merge_bed('outfile_name', ['/bad/path/to/unexisting/file'])

    def test_merge_single_file(self):
        """
        Test merging of single file
        """
        bed1 = [
            ['1', '4', '5', '.', '1', '+'],
            ['1', '5', '6', '.', '1', '+'],
            ['1', '5', '6', '.', '1', '-'],
            ['2', '5', '6', '.', '1', '+'],
            ['1', '5', '6', '.', '2', '+'],
            ['1', '8', '9', '.', '1', '+'],
        ]
        expected = [
            ['1', '4', '5', '.', '1', '+'],
            ['1', '5', '6', '.', '3', '+'],
            ['1', '5', '6', '.', '1', '-'],
            ['1', '8', '9', '.', '1', '+'],
            ['2', '5', '6', '.', '1', '+'],
        ]
        out = merge_bed_wrapper([bed1])

        self.assertEqual(out, expected)

    def test_merge_more_files(self):
        """
        Test that merging in sone well.
        """
        bed1 = [
            ['1', '4', '5', '.', '1', '+'],
            ['1', '5', '6', '.', '1', '+'],
            ['1', '5', '6', '.', '1', '-'],
        ]
        bed2 = [
            ['1', '5', '6', '.', '2', '+'],
            ['2', '5', '6', '.', '1', '+'],
            ['1', '8', '9', '.', '1', '+'],
        ]
        expected = [
            ['1', '4', '5', '.', '1', '+'],
            ['1', '5', '6', '.', '3', '+'],
            ['1', '5', '6', '.', '1', '-'],
            ['1', '8', '9', '.', '1', '+'],
            ['2', '5', '6', '.', '1', '+'],
        ]
        out = merge_bed_wrapper([bed1, bed2])

        self.assertEqual(out, expected)


if __name__ == '__main__':
    unittest.main()
