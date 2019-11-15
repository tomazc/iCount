# pylint: disable=missing-docstring, protected-access
import os
import warnings
import unittest
from unittest.mock import patch  # pylint: disable=unused-import

from pybedtools import create_interval_from_list, BedTool

from iCount.genomes import region
from iCount.genomes.constants import SUBTYPE_GROUPS
from iCount.tests.utils import make_file_from_list, make_list_from_file, get_temp_file_name, get_temp_dir


class TestConstructBorders(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", ResourceWarning)

    def test_basic(self):
        segmentation = [
            # Transcript #1
            ['1', '.', 'ncRNA', '1', '10', '.', '+', '.', 'biotype "A"; gene_name "X";'],
            ['1', '.', 'intron', '11', '20', '.', '+', '.', 'biotype "A"; gene_name "X";'],
            ['1', '.', 'CDS', '21', '30', '.', '+', '.', 'biotype "A"; gene_name "X";'],
            ['1', '.', 'UTR3', '31', '40', '.', '+', '.', 'biotype "A"; gene_name "X";'],
            # Transcript #1
            ['1', '.', 'CDS', '5', '14', '.', '+', '.', 'biotype "A"; gene_name "X";'],
            ['1', '.', 'intron', '15', '24', '.', '+', '.', 'biotype "A"; gene_name "X";'],
            ['1', '.', 'CDS', '25', '34', '.', '+', '.', 'biotype "A"; gene_name "X";'],
            # Also negative strand:
            ['1', '.', 'CDS', '3', '32', '.', '-', '.', 'biotype "A"; gene_name "X";'],
        ]
        expected = [
            ['1', '0', '4', '.', '.', '+'],
            ['1', '4', '10', '.', '.', '+'],
            ['1', '10', '14', '.', '.', '+'],
            ['1', '14', '20', '.', '.', '+'],
            ['1', '20', '24', '.', '.', '+'],
            ['1', '24', '30', '.', '.', '+'],
            ['1', '30', '34', '.', '.', '+'],
            ['1', '34', '40', '.', '.', '+'],
            ['1', '2', '32', '.', '.', '-'],
        ]

        segmentation_file = make_file_from_list(segmentation)
        borders_file = region.construct_borders(BedTool(segmentation_file))
        results = make_list_from_file(borders_file, fields_separator='\t')
        self.assertEqual(
            expected,
            # Sort results by chrom, strand, start, stop
            sorted(results, key=lambda x: (x[0], x[-1], int(x[1]), int(x[2])))
        )


class TestSimplifyBiotype(unittest.TestCase):

    def test_simplify(self):
        self.assertEqual('mRNA', region.simplify_biotype('CDS', 'IG_C_gene'))
        self.assertEqual('pre-mRNA', region.simplify_biotype('intron', 'IG_C_gene'))
        self.assertEqual('lncRNA', region.simplify_biotype('UTR3', 'TEC'))
        self.assertEqual('lncRNA', region.simplify_biotype('ncRNA', 'protein_coding'))

    def test_uniqness_of_entries(self):
        """
        Ensure that entries in SUBTYPE_GROUPS do not repeat.
        """
        all_elements = []
        for _, group_elements in SUBTYPE_GROUPS.items():
            all_elements.extend(group_elements)

        self.assertEqual(len(all_elements), len(set(all_elements)))


class TestMakeUniqRegion(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", ResourceWarning)

    def test_basic(self):
        # seg is composed of borders(BED6) and segment(GTF) interval:
        seg = create_interval_from_list(
            ['1', '0', '10', '.', '.', '+'] + ['.', '.', '.', '.', '.', '.', '.', '.', '.'])
        types = ['UTR3']
        subtypes = ['TEC']
        genes = [('id1', 'A', 50)]

        interval = region.make_uniq_region(seg, types, subtypes, genes)
        self.assertEqual(interval[:], ['1', '.', 'UTR3', '1', '10', '.', '+', '.',
                                       'gene_id "id1"; gene_name "A"; biotype "lncRNA";'])

    def test_highest_rated_type(self):
        # seg is compositon of BED6 and GTF interval:
        seg = create_interval_from_list(
            ['1', '0', '10', '.', '.', '+'] + ['.', '.', '.', '.', '.', '.', '.', '.', '.'])
        types = ['UTR3', 'intron', 'UTR5']
        subtypes = ['protein_coding', 'TEC', 'non_stop_decay']
        genes = [('id1', 'A', 20), ('id1', 'A', 20), ('id1', 'A', 20)]

        interval = region.make_uniq_region(seg, types, subtypes, genes)
        self.assertEqual(interval[:], ['1', '.', 'UTR3', '1', '10', '.', '+', '.',
                                       'gene_id "id1"; gene_name "A"; biotype "mRNA";'])

    def test_multiple_biotypes(self):
        # seg is compositon of BED6 and GTF interval:
        seg = create_interval_from_list(
            ['1', '0', '10', '.', '.', '+'] + ['.', '.', '.', '.', '.', '.', '.', '.', '.'])
        types = ['intron', 'intron']
        subtypes = ['protein_coding', 'TEC']
        genes = [('id1', 'A', 20), ('id1', 'A', 20)]

        interval = region.make_uniq_region(seg, types, subtypes, genes)
        self.assertEqual(interval[:], ['1', '.', 'intron', '1', '10', '.', '+', '.',
                                       'gene_id "id1"; gene_name "A"; biotype "lncRNA,pre-mRNA";'])

    def test_take_longer_gene(self):
        # seg is compositon of BED6 and GTF interval:
        seg = create_interval_from_list(
            ['1', '0', '10', '.', '.', '+'] + ['.', '.', '.', '.', '.', '.', '.', '.', '.'])
        types = ['CDS', 'CDS']
        subtypes = ['protein_coding', 'protein_coding']
        genes = [('id1', 'A', 20), ('id2', 'B', 40)]

        interval = region.make_uniq_region(seg, types, subtypes, genes)
        self.assertEqual(interval[:], ['1', '.', 'CDS', '1', '10', '.', '+', '.',
                                       'gene_id "id2"; gene_name "B"; biotype "mRNA";'])

    def test_utr3(self):
        # seg is compositon of BED6 and GTF interval:
        seg = create_interval_from_list(
            ['1', '0', '10', '.', '.', '+'] + ['.', '.', '.', '.', '.', '.', '.', '.', '.'])
        types = ['intron']
        subtypes = ['3prime_overlapping_ncRNA']
        genes = [('id1', 'A', 20)]

        interval = region.make_uniq_region(seg, types, subtypes, genes)
        self.assertEqual(interval[:], ['1', '.', 'UTR3', '1', '10', '.', '+', '.',
                                       'gene_id "id1"; gene_name "A"; biotype "mRNA";'])

    def test_intergenic(self):
        # seg is compositon of BED6 and GTF interval:
        seg = create_interval_from_list(
            ['1', '0', '10', '.', '.', '+'] + ['.', '.', '.', '.', '.', '.', '.', '.', '.'])
        types = ['intergenic']
        subtypes = [None]
        genes = [('.', None, 0)]

        interval = region.make_uniq_region(seg, types, subtypes, genes)
        self.assertEqual(interval[:], ['1', '.', 'intergenic', '1', '10', '.', '+', '.',
                                       'gene_id "."; gene_name "None"; biotype "";'])


class TestMergeRegions(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", ResourceWarning)
        self.tmp = get_temp_file_name()

    def test_basic(self):
        # seg is compositon of BED6 and GTF interval:
        nonmerged = make_file_from_list([
            ['1', '.', 'UTR3', '1', '10', '.', '+', '.', 'biotype "lncRNA";gene_id "id1";'],
            ['1', '.', 'UTR3', '11', '20', '.', '+', '.', 'biotype "lncRNA";gene_id "id1";'],
            ['1', '.', 'UTR3', '21', '30', '.', '+', '.', 'biotype "lncRNA";gene_id "id2";'],
            ['1', '.', 'UTR3', '31', '40', '.', '+', '.', 'biotype "lncRNA";gene_id "id1";'],
            ['1', '.', 'UTR3', '31', '40', '.', '-', '.', 'biotype "lncRNA";gene_id "id1";'],
        ])

        expected = [
            ['1', '.', 'UTR3', '1', '20', '.', '+', '.', 'biotype "lncRNA";gene_id "id1";'],
            ['1', '.', 'UTR3', '21', '30', '.', '+', '.', 'biotype "lncRNA";gene_id "id2";'],
            ['1', '.', 'UTR3', '31', '40', '.', '+', '.', 'biotype "lncRNA";gene_id "id1";'],
            ['1', '.', 'UTR3', '31', '40', '.', '-', '.', 'biotype "lncRNA";gene_id "id1";'],
        ]

        region.merge_regions(nonmerged, self.tmp)
        results = make_list_from_file(self.tmp, fields_separator='\t')
        # Since order of attrs can be arbitrary, equality checks are more complex:
        for res, exp in zip(results, expected):
            self.assertEqual(res[:8], exp[:8])
            self.assertEqual(
                ';'.join(sorted(res[8].split(';'))),
                ';'.join(sorted(exp[8].split(';'))),
            )


class TestSummaryTemplates(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", ResourceWarning)

    def test_templates1(self):
        out_dir = get_temp_dir()
        segmentation = make_file_from_list([
            ['1', '.', 'intergenic', '1', '10', '.', '+', '.', 'gene_id ".";'],
            ['1', '.', 'UTR3', '11', '20', '.', '+', '.', 'biotype "mRNA";gene_name "ABC";gene_id "G1";'],
            ['1', '.', 'intron', '21', '30', '.', '+', '.', 'biotype "lncRNA";gene_name "ABC";gene_id "G1";'],
            ['1', '.', 'CDS', '31', '40', '.', '+', '.', 'biotype "mRNA";gene_name "DEF";gene_id "G2";'],
            ['1', '.', 'intron', '41', '50', '.', '+', '.', 'biotype "sRNA,lncRNA";gene_name "DEF"; gene_id "G2";'],
        ])
        region.summary_templates(segmentation, out_dir)

        results_type = make_list_from_file(os.path.join(out_dir, region.TEMPLATE_TYPE), '\t')
        self.assertEqual(results_type, [
            ['CDS', '10'],
            ['UTR3', '10'],
            ['intron', '20'],
            ['intergenic', '10'],
        ])

        results_subtype = make_list_from_file(os.path.join(out_dir, region.TEMPLATE_SUBTYPE), fields_separator='\t')
        self.assertEqual(results_subtype, [
            ['CDS mRNA', '10'],
            ['UTR3 mRNA', '10'],
            ['intron lncRNA', '15'],
            ['intron sRNA', '5'],
            ['intergenic', '10'],
        ])

        results_gene = make_list_from_file(os.path.join(out_dir, region.TEMPLATE_GENE), fields_separator='\t')
        self.assertEqual(results_gene, [
            ['.', '', '10'],
            ['G1', 'ABC', '20'],
            ['G2', 'DEF', '20'],
        ])


class TestMakeRegionsFile(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", ResourceWarning)
        self.dir = get_temp_dir()

    def test_basic(self):
        segmentation = [
            ['1', '.', 'gene', '1', '50', '.', '+', '.', 'biotype "lincRNA"; gene_name "A"; gene_id "X";'],
            # Transcript #1
            ['1', '.', 'transcript', '1', '40', '.', '+', '.', 'biotype "lincRNA"; gene_name "A"; gene_id "X";'],
            ['1', '.', 'ncRNA', '1', '10', '.', '+', '.', 'biotype "lincRNA"; gene_name "A"; gene_id "X";'],
            ['1', '.', 'UTR5', '11', '20', '.', '+', '.', 'biotype "lincRNA"; gene_name "A"; gene_id "X";'],
            ['1', '.', 'CDS', '21', '30', '.', '+', '.', 'biotype "lincRNA"; gene_name "A"; gene_id "X";'],
            ['1', '.', 'intron', '31', '35', '.', '+', '.', 'biotype "lincRNA"; gene_name "A"; gene_id "X";'],
            ['1', '.', 'CDS', '36', '40', '.', '+', '.', 'biotype "lincRNA"; gene_name "A"; gene_id "X";'],
            # Transcript #2
            ['1', '.', 'transcript', '10', '50', '.', '+', '.', 'biotype "rRNA"; gene_name "A"; gene_id "X";'],
            ['1', '.', 'ncRNA', '10', '18', '.', '+', '.', 'biotype "rRNA"; gene_name "A"; gene_id "X";'],
            ['1', '.', 'UTR5', '19', '25', '.', '+', '.', 'biotype "rRNA"; gene_name "A"; gene_id "X";'],
            ['1', '.', 'CDS', '26', '32', '.', '+', '.', 'biotype "rRNA"; gene_name "A"; gene_id "X";'],
            ['1', '.', 'intron', '33', '39', '.', '+', '.', 'biotype "rRNA"; gene_name "A"; gene_id "X";'],
            ['1', '.', 'CDS', '40', '44', '.', '+', '.', 'biotype "rRNA"; gene_name "A"; gene_id "X";'],
            ['1', '.', 'UTR3', '45', '50', '.', '+', '.', 'biotype "rRNA"; gene_name "A"; gene_id "X";'],
            # Itergenic
            ['1', '.', 'intergenic', '51', '100', '.', '+', '.', 'gene_id ".";'],

        ]
        expected = [
            ['1', '.', 'ncRNA', '1', '9', '.', '+', '.', 'gene_id "X";biotype "lncRNA";gene_name "A";'],
            ['1', '.', 'ncRNA', '10', '10', '.', '+', '.', 'gene_id "X";biotype "lncRNA,rRNA";gene_name "A";'],
            ['1', '.', 'UTR5', '11', '18', '.', '+', '.', 'gene_id "X";biotype "lncRNA";gene_name "A";'],
            ['1', '.', 'UTR5', '19', '20', '.', '+', '.', 'gene_id "X";biotype "lncRNA,rRNA";gene_name "A";'],
            ['1', '.', 'CDS', '21', '25', '.', '+', '.', 'gene_id "X";biotype "lncRNA";gene_name "A";'],
            ['1', '.', 'CDS', '26', '30', '.', '+', '.', 'gene_id "X";biotype "lncRNA,rRNA";gene_name "A";'],
            ['1', '.', 'CDS', '31', '32', '.', '+', '.', 'gene_id "X";biotype "rRNA";gene_name "A";'],
            ['1', '.', 'intron', '33', '35', '.', '+', '.', 'gene_id "X";biotype "lncRNA,rRNA";gene_name "A";'],
            ['1', '.', 'CDS', '36', '39', '.', '+', '.', 'gene_id "X";biotype "lncRNA";gene_name "A";'],
            ['1', '.', 'CDS', '40', '40', '.', '+', '.', 'gene_id "X";biotype "lncRNA,rRNA";gene_name "A";'],
            ['1', '.', 'CDS', '41', '44', '.', '+', '.', 'gene_id "X";biotype "rRNA";gene_name "A";'],
            ['1', '.', 'UTR3', '45', '50', '.', '+', '.', 'gene_id "X";biotype "rRNA";gene_name "A";'],
            ['1', '.', 'intergenic', '51', '100', '.', '+', '.', 'gene_id ".";biotype "";gene_name "None";'],
        ]

        segmentation_file = make_file_from_list(segmentation, sort=True)
        region.make_regions(segmentation_file, self.dir)
        results = make_list_from_file(os.path.join(self.dir, region.REGIONS_FILE), fields_separator='\t')

        # Since order of attrs can be arbitrary, equality checks are more complex:
        for res, exp in zip(results, expected):
            self.assertEqual(res[:8], exp[:8])
            self.assertEqual(
                ';'.join(sorted(res[8].split(';'))),
                ';'.join(sorted(exp[8].split(';'))),
            )


if __name__ == '__main__':
    unittest.main()
