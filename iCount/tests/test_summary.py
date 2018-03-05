# pylint: disable=missing-docstring, protected-access
import os
import unittest
import warnings

from iCount.analysis import summary
from iCount.genomes import segment
from iCount.tests.utils import make_file_from_list, make_list_from_file, get_temp_dir


class TestMakeSummaryReport(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", ResourceWarning)
        self.out_dir = get_temp_dir()
        self.type_header = ['Type', 'Length', 'cDNA #', 'cDNA %']
        self.subtype_header = ['Subtype', 'Length', 'cDNA #', 'cDNA %']
        self.gene_header = ['Gene name (Gene ID)', 'Length', 'cDNA #', 'cDNA %']

    def get_summary_reports(self, annotation, cross_links):
        """Help running tests for ``summary_report`` with less clutter."""
        annotation_file = make_file_from_list(annotation)
        cross_links_file = make_file_from_list(cross_links)

        segment.summary_templates(annotation_file, self.out_dir)
        summary.summary_reports(annotation_file, cross_links_file, self.out_dir, self.out_dir)
        return [
            make_list_from_file(os.path.join(self.out_dir, segment.SUMMARY_TYPE), '\t'),
            make_list_from_file(os.path.join(self.out_dir, segment.SUMMARY_SUBTYPE), '\t'),
            make_list_from_file(os.path.join(self.out_dir, segment.SUMMARY_GENE), '\t'),
        ]

    def test_diff_chromosome_naming(self):
        """
        Exception is raised if chromosome naming is inconsistent.
        """
        cross_links = [
            ['chr1', '15', '16', '.', '5', '+'],
        ]
        annotation = [
            ['1', '.', 'CDS', '1', '10', '.', '+', '.', 'biotype "mRNA";gene_id "G1";'],
        ]
        message = r"No intersections found. This may be caused by .*"
        with self.assertRaisesRegex(ValueError, message):
            self.get_summary_reports(annotation, cross_links)

    def test_diff_only_strand1(self):
        """
        Same coords but diff strand and same type.
        """
        cross_links = [
            ['1', '5', '6', '.', '3', '+'],
            ['1', '5', '6', '.', '2', '-']]
        annotation = [
            ['1', '.', 'CDS', '1', '10', '.', '+', '.', 'biotype "lncRNA";gene_id "G1";gene_name "ABC";'],
            ['1', '.', 'CDS', '1', '10', '.', '-', '.', 'biotype "lncRNA";gene_id "G1";gene_name "ABC";'],
        ]

        results_type, results_subtype, results_gene = self.get_summary_reports(annotation, cross_links)

        self.assertEqual(results_type, [
            self.type_header,
            ['CDS', '20', '5', '100.0'],
        ])

        self.assertEqual(results_subtype, [
            self.subtype_header,
            ['CDS lncRNA', '20', '5', '100.0'],
        ])

        self.assertEqual(results_gene, [
            self.gene_header,
            ['ABC (G1)', '20', '5', '100.0'],
        ])

    def test_diff_only_strand2(self):
        """
        Same coords but diff strand and diff type.
        """
        cross_links = [
            ['1', '5', '6', '.', '3', '+'],
            ['1', '5', '6', '.', '1', '-']]
        annotation = [
            ['1', '.', 'CDS', '1', '10', '.', '+', '.', 'biotype "lncRNA";gene_id ".";'],
            ['1', '.', 'CDS', '1', '10', '.', '-', '.', 'biotype "miRNA";gene_id ".";'],
        ]

        results_type, results_subtype, results_gene = self.get_summary_reports(annotation, cross_links)

        self.assertEqual(results_type, [
            self.type_header,
            ['CDS', '20', '4', '100.0'],
        ])

        self.assertEqual(results_subtype, [
            self.subtype_header,
            ['CDS lncRNA', '10', '3', '75.0'],
            ['CDS miRNA', '10', '1', '25.0'],
        ])

        self.assertEqual(results_gene, [
            self.gene_header,
            ['intergenic (.)', '20', '4', '100.0'],
        ])

    def test_many_regions(self):
        """
        Multiple annotation regions intersecting one crosslink.
        """
        cross_links = [
            ['1', '5', '6', '.', '2', '+'],
            ['1', '15', '16', '.', '2', '+'],
        ]
        annotation = [
            ['1', '.', 'CDS', '1', '10', '.', '+', '.', 'biotype "rRNA,sRNA";gene_id ".";'],
            ['1', '.', 'ncRNA', '11', '30', '.', '+', '.', 'biotype "rRNA,snRNA";gene_id ".";'],
        ]

        results_type, results_subtype, _ = self.get_summary_reports(annotation, cross_links)

        self.assertEqual(results_type, [
            self.type_header,
            ['CDS', '10', '2', '50.0'],
            ['ncRNA', '20', '2', '50.0'],
        ])

        self.assertEqual(results_subtype, [
            self.subtype_header,
            ['CDS rRNA', '5', '1', '25.0'],
            ['CDS sRNA', '5', '1', '25.0'],
            ['ncRNA rRNA', '10', '1', '25.0'],
            ['ncRNA snRNA', '10', '1', '25.0'],
        ])

    def test_any_annotation(self):
        """
        Multiple annotation regions intersecting one crosslink.
        """
        cross_links = [
            ['1', '5', '6', '.', '2', '+'],
            ['1', '15', '16', '.', '2', '+'],
        ]
        annotation = [
            ['1', '.', 'XYZ', '1', '10', '.', '+', '.', 'biotype "magic";gene_name "ABC";gene_id "123";'],
            ['1', '.', 'WHY', '11', '30', '.', '+', '.', 'biotype "poison";gene_name "DEF";gene_id "456";'],
        ]

        results_type, results_subtype, results_gene = self.get_summary_reports(annotation, cross_links)

        self.assertEqual(sorted(results_type), sorted([
            self.type_header,
            ['WHY', '20', '2', '50.0'],
            ['XYZ', '10', '2', '50.0'],
        ]))

        self.assertEqual(sorted(results_subtype), sorted([
            self.subtype_header,
            ['WHY poison', '20', '2', '50.0'],
            ['XYZ magic', '10', '2', '50.0'],
        ]))

        self.assertEqual(sorted(results_gene), sorted([
            self.gene_header,
            ['ABC (123)', '10', '2', '50.0'],
            ['DEF (456)', '20', '2', '50.0'],
        ]))


if __name__ == '__main__':
    unittest.main()
