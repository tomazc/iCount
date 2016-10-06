import unittest
import warnings

import iCount
from iCount.tests.utils import make_file_from_list, make_list_from_file, get_temp_file_name


class TestGetGenes(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", ResourceWarning)

    def test_all_good(self):

        gtf = make_file_from_list([
            ['1', '.', 'transcript', '200', '250', '.', '+', '.',
             'gene_id "G1"; transcript_id "T1";'],
            ['1', '.', 'transcript', '100', '300', '.', '+', '.',
             'gene_id "G1"; transcript_id "T2";'],
            ['2', '.', 'gene', '400', '500', '.', '+', '.', 'gene_id "G3";'],
        ])

        result = make_list_from_file(
            iCount.genomes.genes.get_genes(gtf, get_temp_file_name()), fields_separator='\t')

        expected = [
            ['1', '.', 'gene', '99', '300', '.', '+', '.', 'gene_id "G1"; transcript_id "T1";'],
            ['2', '.', 'gene', '399', '500', '.', '+', '.', 'gene_id "G3";'],
        ]

        self.assertEqual(result, expected)
