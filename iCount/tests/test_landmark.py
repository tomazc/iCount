# pylint: disable=missing-docstring, protected-access
import unittest
import warnings

from pybedtools import create_interval_from_list

from iCount.genomes import landmark
from iCount.tests.utils import get_temp_file_name, make_file_from_list, make_list_from_file


class TestGetGeneName(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", (ResourceWarning, ImportWarning))

    def test_basic(self):
        seg_level0 = create_interval_from_list(['1', '.', 'CDS', '1', '2', '.', '+', '.', 'gene_name "G0";'])
        seg_level1 = create_interval_from_list(['1', '.', 'UTR3', '1', '2', '.', '+', '.', 'gene_name "G1";'])
        seg_level4 = create_interval_from_list(['1', '.', 'intron', '1', '2', '.', '+', '.', 'gene_name "G2";'])

        self.assertEqual('G0', landmark.get_gene_name(seg_level0, seg_level1))
        self.assertEqual('G0', landmark.get_gene_name(seg_level1, seg_level0))
        self.assertEqual('G1', landmark.get_gene_name(seg_level1, seg_level4))
        with self.assertRaises(ValueError):
            self.assertEqual('B', landmark.get_gene_name(seg_level0, seg_level0))


class TestMakeSingleTypeLandmarks(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", (ResourceWarning, ImportWarning))

    def test_basic(self):
        regions = make_file_from_list([
            ['chr1', '.', 'CDS', '150', '200', '.', '+', '.', 'gene_name "A";'],
            ['chr1', '.', 'intron', '201', '351', '.', '+', '.', 'gene_name "A";'],
        ])
        landmarks = landmark.make_single_type_landmarks(regions, 'exon-intron')
        landmarks = [list(map(str, item)) for item in landmarks]
        self.assertEqual(landmarks, [
            ['chr1', '200', '201', 'exon-intron;A', '.', '+'],
        ])

        regions = make_file_from_list([
            ['chr1', '.', 'CDS', '150', '200', '.', '-', '.', 'gene_name "A";'],
            ['chr1', '.', 'intron', '201', '351', '.', '-', '.', 'gene_name "A";'],
        ])
        landmarks = landmark.make_single_type_landmarks(regions, 'intron-exon')
        landmarks = [list(map(str, item)) for item in landmarks]
        self.assertEqual(landmarks, [
            ['chr1', '199', '200', 'intron-exon;A', '.', '-'],
        ])

    def test_limits_upstream(self):
        """Landmarks with too short upstream segment should not be used."""
        regions = make_file_from_list([
            ['chr1', '.', 'CDS', '151', '200', '.', '+', '.', 'gene_name "A";'],
            ['chr1', '.', 'intron', '201', '400', '.', '+', '.', 'gene_name "A";'],
        ])
        landmarks = landmark.make_single_type_landmarks(regions, 'exon-intron')
        landmarks = [list(map(str, item)) for item in landmarks]
        self.assertEqual(landmarks, [])

        regions = make_file_from_list([
            ['chr1', '.', 'CDS', '150', '200', '.', '-', '.', 'gene_name "A";'],
            ['chr1', '.', 'intron', '201', '350', '.', '-', '.', 'gene_name "A";'],
        ])
        landmarks = landmark.make_single_type_landmarks(regions, 'intron-exon')
        landmarks = [list(map(str, item)) for item in landmarks]
        self.assertEqual(landmarks, [])

    def test_limits_downstream(self):
        """Landmarks with too short upstream segment should not be used."""
        regions = make_file_from_list([
            ['chr1', '.', 'CDS', '150', '200', '.', '+', '.', 'gene_name "A";'],
            ['chr1', '.', 'intron', '201', '350', '.', '+', '.', 'gene_name "A";'],
        ])
        landmarks = landmark.make_single_type_landmarks(regions, 'exon-intron')
        landmarks = [list(map(str, item)) for item in landmarks]
        self.assertEqual(landmarks, [])

        regions = make_file_from_list([
            ['chr1', '.', 'CDS', '151', '200', '.', '-', '.', 'gene_name "A";'],
            ['chr1', '.', 'intron', '201', '351', '.', '-', '.', 'gene_name "A";'],
        ])
        landmarks = landmark.make_single_type_landmarks(regions, 'intron-exon')
        landmarks = [list(map(str, item)) for item in landmarks]
        self.assertEqual(landmarks, [])


class TestMakeLandmarks(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", (ResourceWarning, ImportWarning))

    def test_basic(self):
        regions = make_file_from_list([
            ['chr1', '.', 'CDS', '150', '200', '.', '+', '.', 'gene_name "A";'],
            ['chr1', '.', 'intron', '201', '400', '.', '+', '.', 'gene_name "A";'],
            ['chr1', '.', 'CDS', '401', '600', '.', '+', '.', 'gene_name "A";'],
        ])

        landmarks = get_temp_file_name(extension='bed')
        landmark.make_landmarks(regions, landmarks)
        self.assertEqual(make_list_from_file(landmarks), [
            ['chr1', '200', '201', 'exon-intron;A', '.', '+'],
            ['chr1', '400', '401', 'intron-exon;A', '.', '+'],
        ])


if __name__ == '__main__':
    unittest.main()
