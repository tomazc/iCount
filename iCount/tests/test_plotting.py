# pylint: disable=missing-docstring, protected-access
import os
import unittest
import warnings

from iCount.plotting import rnamap, rnaheatmap, rnacombined
from iCount.tests.utils import get_temp_file_name, make_file_from_list


def get_example_results_file():
    """Produce some thest file."""
    def zeros(size):
        """Return array of zeros of size ``size``."""
        return [0 for i in range(size)]

    fname = make_file_from_list(bedtool=False, data=[
        ['total_cdna:1000000'],
        ['.'] + [str(pos) for pos in range(-50, 51)],
        ['chr1_+_100'] + zeros(51) + ['3', '5'] + zeros(48),
        ['chr1_+_200'] + zeros(49) + ['3', '0', '0', '1'] + zeros(48),
    ])

    return fname


class TestPlotRnaMap(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", (ResourceWarning, ImportWarning))

        self.results_file = get_example_results_file()

    def test_normalize_cpm(self):

        normalized = rnamap.normalize_cpm(1, 10**6)
        self.assertEqual(normalized, 1)

        normalized = rnamap.normalize_cpm(2, 10**6 * 2)
        self.assertEqual(normalized, 1)

    def test_parse(self):
        data, landmark_count = rnamap.parse_results(self.results_file)
        self.assertEqual(landmark_count, 2)

        expected = {
            -1: 3.0,
            1: 3.0,
            2: 6.0,
        }
        for pos, score in data.items():
            self.assertEqual(score, expected.get(pos, 0.0))

    def test_plot(self):
        outfile = get_temp_file_name(extension='png')
        rnamap.plot_rnamap(self.results_file, outfile)
        self.assertTrue(os.path.isfile(outfile))


class TestSmooth(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", (ResourceWarning, ImportWarning))

    def test_smoothe(self):
        # pylint: disable=missing-docstring, protected-access
        values = [0, 2, 0]
        self.assertEqual(rnamap.smooth(values, 1), [1, 2 / 3, 1])

        values = [0, 2, 4, 1, 0]
        self.assertEqual(rnamap.smooth(values, 1), [1, 2, 7 / 3, 5 / 3, 0.5])
        self.assertEqual(rnamap.smooth(values, 2), [2, 7 / 4, 7 / 5, 7 / 4, 5 / 3])


class TestPlotRnaHeatMap(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", (ResourceWarning, ImportWarning))

        self.results_file = get_example_results_file()

    def test_make_position_to_bin1(self):
        pos2bins = rnaheatmap.make_position_to_bin([0, 2, 4])
        self.assertEqual(pos2bins, {0: 0, 1: 0, 2: 2, 3: 2, 4: 2})

    def test_parse(self):
        data = rnaheatmap.parse_results(self.results_file, -50, 50, top_n=2, binsize=10)
        self.assertEqual(list(data.index), ['chr1_+_100', 'chr1_+_200'])
        self.assertEqual(list(data.columns), list(range(-50, 50, 10)))
        self.assertEqual(data.at['chr1_+_100', 0], 8.0)
        self.assertEqual(data.at['chr1_+_200', -10], 3.0)
        self.assertEqual(data.at['chr1_+_200', 0], 1.0)

    def test_plot(self):
        outfile = get_temp_file_name(extension='png')
        rnaheatmap.plot_rnaheatmap(self.results_file, outfile, top_n=2, binsize=10)
        self.assertTrue(os.path.isfile(outfile))


class TestPlotRnaCombined(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", (ResourceWarning, ImportWarning))

        self.results_file = get_example_results_file()

    def test_plot(self):
        outfile = get_temp_file_name(extension='png')
        rnacombined.plot_combined(self.results_file, outfile, top_n=2, nbins=50)
        self.assertTrue(os.path.isfile(outfile))
