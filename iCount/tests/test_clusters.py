import unittest

from pybedtools import create_interval_from_list

from iCount.analysis import clusters
from iCount.tests.utils import make_file_from_list, make_list_from_file, \
    get_temp_file_name


class TestClusters(unittest.TestCase):

    def test_fix_proper_bed6_format(self):
        feature = create_interval_from_list(['1', '1', '10', '+', '5'])

        converted = clusters._fix_proper_bed6_format(feature)
        expected = create_interval_from_list(['1', '1', '10', '', '5', '+'])
        self.assertEqual(expected, converted)

    def test_clusters(self):
        fin_sites = make_file_from_list([
            ['1', '4', '5', '.', '1', '+'],
            ['1', '4', '5', '.', '1', '-'],
            ['1', '5', '6', '.', '1', '+'],
            ['1', '10', '11', '.', '1', '+'],
        ])

        fout_clusters = get_temp_file_name()

        result = make_list_from_file(clusters.run(
            fin_sites, fout_clusters, dist=3), fields_separator='\t')

        expected = [
            ['1', '4', '6', '.', '2', '+'],
            ['1', '4', '5', '.', '1', '-'],
            ['1', '10', '11', '.', '1', '+'],
        ]

        self.assertEqual(expected, result)


if __name__ == '__main__':
    unittest.main()
