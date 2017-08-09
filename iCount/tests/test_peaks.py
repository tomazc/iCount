# pylint: disable=missing-docstring, protected-access

import unittest
import warnings

from iCount.analysis import peaks
from iCount.tests.utils import get_temp_file_name, make_file_from_list, \
    make_list_from_file


class TestPeaks(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", ResourceWarning)

    def test_sum_within_window(self):
        sites1 = [
            (10, 1), (11, 1), (12, 1), (20, 2),
        ]
        summed_sites_w1 = [
            (10, 2), (11, 3), (12, 2), (20, 2),
        ]
        sites2 = [
            (10, 1), (11, 1), (12, 1), (13, 1), (14, 1), (15, 1), (20, 2),
        ]
        summed_sites_w3 = [
            (10, 4), (11, 5), (12, 6), (13, 6), (14, 5), (15, 4), (20, 2),
        ]

        self.assertEqual(peaks._sum_within_window([]), [])

        self.assertEqual(
            peaks._sum_within_window(sites1, half_window=1), summed_sites_w1)
        self.assertEqual(
            peaks._sum_within_window(sites2, half_window=3), summed_sites_w3)

        # shuffle the input:
        sites3 = [
            (12, 1), (10, 1), (20, 2), (11, 1),
        ]
        summed_sites_w1_2 = [
            (12, 2), (10, 2), (20, 2), (11, 3),
        ]
        self.assertEqual(
            peaks._sum_within_window(sites3, half_window=1), summed_sites_w1_2)

    def test_sum_within_window_nopos(self):
        sites = [
            (10, 1), (11, 1), (12, 1), (13, 1), (14, 1), (15, 1), (20, 2),
        ]
        summed_sites_w1 = [2, 3, 3, 3, 3, 2, 2]
        summed_sites_w3 = [4, 5, 6, 6, 5, 4, 2]

        self.assertEqual(peaks._sum_within_window_nopos([]), [])

        self.assertEqual(
            peaks._sum_within_window_nopos(sites, half_window=1), summed_sites_w1)
        self.assertEqual(
            peaks._sum_within_window_nopos(sites, half_window=3), summed_sites_w3)

    def test_cumulative_prob(self):
        vals = [2, 3, 3, 3, 2]
        max_val = 5

        expected = [1., 1., 1., 0.6, 0., 0.]
        result = list(peaks.cumulative_prob(vals, max_val))
        self.assertEqual(expected, result)

    def test_get_avg_rnd_distrib(self):
        size = 5
        total_hits = 5
        perms = 10000

        expected = [1., 1., 1., 0.86, 0.63, 0.30]
        result = peaks.get_avg_rnd_distrib(size, total_hits, 1, perms=perms)
        # Repeat the same call, so also lines that use cache are executed:
        result = peaks.get_avg_rnd_distrib(size, total_hits, 1, perms=perms)

        for res, exp, in zip(result, expected):
            self.assertAlmostEqual(res, exp, delta=0.02)

    def test_run(self):
        fin_annotation = make_file_from_list([
            ['1', '.', 'gene', '10', '20', '.', '+', '.', 'gene_name "A"; gene_id "1";'],
            ['1', '.', 'transcript', '10', '20', '.', '+', '.', 'gene_name "B"; gene_id "1";'],
            ['2', '.', 'CDS', '10', '20', '.', '+', '.', 'gene_name "C"; gene_id "1";'],
        ])

        fin_sites = make_file_from_list([
            ['1', '14', '15', '.', '3', '+'],
            ['1', '16', '17', '.', '5', '+'],
            ['2', '16', '17', '.', '5', '+'],
        ])

        fout_peaks = get_temp_file_name(extension='.bed.gz')
        fout_scores = get_temp_file_name(extension='.tsv.gz')

        peaks.run(fin_annotation, fin_sites, fout_peaks,
                  scores=fout_scores)

        out_peaks = make_list_from_file(fout_peaks, fields_separator='\t')
        out_scores = make_list_from_file(fout_scores, fields_separator='\t')
        # Remove header:
        out_scores = out_scores[1:]

        expected_peaks = [
            ['1', '14', '15', 'A-1', '3', '+'],
            ['1', '16', '17', 'A-1', '5', '+'],
        ]
        expected_scores = [
            ['1', '14', '+', 'A', '1', '3', '8', '0.036198'],
            ['1', '16', '+', 'A', '1', '5', '8', '0.036198'],
            ['2', '16', '+', 'not_annotated', 'not_annotated', '5', 'not_calculated', '1'],
        ]

        self.assertEqual(out_peaks, expected_peaks)
        self.assertEqual(out_scores, expected_scores)


if __name__ == '__main__':
    unittest.main()
