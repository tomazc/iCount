# pylint: disable=missing-docstring, protected-access

import unittest
import warnings

from pybedtools import create_interval_from_list

from iCount.analysis import clusters
from iCount.tests.utils import make_file_from_list, make_list_from_file, \
    get_temp_file_name


class TestClusters(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", ResourceWarning)

    def test_fix_proper_bed6_format(self):
        feature = create_interval_from_list(['1', '1', '10', '+', '5'])

        converted = clusters._fix_bed6(feature)
        expected = create_interval_from_list(['1', '1', '10', '', '5', '+'])
        self.assertEqual(expected, converted)

    def test_clusters(self):
        fin_sites = make_file_from_list([
            ['1', '1', '2', '.', '1', '+'],
            ['1', '2', '3', '.', '1', '+'],
            ['1', '3', '4', '.', '1', '+'],
            ['1', '4', '5', '.', '2', '+'],
            ['1', '4', '5', '.', '1', '-'],
            ['1', '5', '6', '.', '1', '+'],
            ['1', '6', '7', '.', '1', '-'],
            ['1', '7', '8', '.', '1', '-'],
            ['1', '10', '11', '.', '1', '+'],
            ['1', '11', '12', '.', '2', '+'],
            ['1', '12', '13', '.', '1', '+'],
        ])

        fin_peaks = make_file_from_list([
            ['1', '4', '5', 'cl1', '1', '+'],
            ['1', '4', '5', 'cl2', '1', '-'],
            ['1', '5', '6', 'cl3', '1', '+'],
            ['1', '11', '12', 'cl4', '2', '+'],
        ])

        fout_clusters = get_temp_file_name()

        clusters.run(fin_sites, fin_peaks, fout_clusters, dist=3, slop=2)
        result = make_list_from_file(fout_clusters, fields_separator='\t')

        expected = [
            ['1', '2', '6', 'cl1,cl3', '5', '+'],
            ['1', '4', '7', 'cl2', '2', '-'],
            ['1', '10', '13', 'cl4', '4', '+'],
        ]

        self.assertEqual(expected, result)


if __name__ == '__main__':
    unittest.main()
