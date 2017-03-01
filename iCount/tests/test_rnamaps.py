# pylint: disable=missing-docstring, protected-access
import os
import unittest
import warnings

from iCount.analysis import rnamaps
from iCount.tests.utils import get_temp_file_name, make_bam_file, make_file_from_list, \
    list_to_intervals, intervals_to_list, make_list_from_file


class TestRun(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", (ResourceWarning, ImportWarning))
        self.gtf_data = list_to_intervals([
            ['1', '.', 'intergenic', '1', '99', '.', '+', '.',
             'gene_id "."; transcript_id ".";'],
            # Gene #1:
            ['1', '.', 'gene', '100', '499', '.', '+', '.',
             'gene_id "G1";'],
            # Transcript #1
            ['1', '.', 'transcript', '100', '249', '.', '+', '.',
             'gene_id "G1"; transcript_id "T1";'],
            ['1', '.', 'UTR5', '100', '149', '.', '+', '.',
             'gene_id "G1"; transcript_id "T1"; exon_number "1";'],
            ['1', '.', 'intron', '150', '199', '.', '+', '.',
             'gene_id "G1"; transcript_id "T1";'],
            ['1', '.', 'CDS', '200', '229', '.', '+', '.',
             'gene_id "G1"; transcript_id "T1"; exon_number "2";'],
            ['1', '.', 'intron', '230', '239', '.', '+', '.',
             'gene_id "G1"; transcript_id "T1";'],
            ['1', '.', 'UTR3', '240', '249', '.', '+', '.',
             'gene_id "G1"; transcript_id "T1"; exon_number "3";'],

            # Transcript #2
            ['1', '.', 'transcript', '240', '499', '.', '+', '.',
             'gene_id "G1"; transcript_id "T2";'],
            ['1', '.', 'CDS', '240', '299', '.', '+', '.',
             'gene_id "G1"; transcript_id "T2"; exon_number "1";'],
            ['1', '.', 'intron', '300', '399', '.', '+', '.',
             'gene_id "G1"; transcript_id "T2";'],
            ['1', '.', 'CDS', '400', '499', '.', '+', '.',
             'gene_id "G1"; transcript_id "T2"; exon_number "2";'],

            # intergenic
            ['1', '.', 'intergenic', '500', '599', '.', '+', '.',
             'gene_id "."; transcript_id ".";'],

            # Gene #1:
            ['1', '.', 'gene', '600', '999', '.', '+', '.',
             'gene_id "G2";'],

            # Transcript #3
            ['1', '.', 'transcript', '600', '799', '.', '+', '.',
             'gene_id "G2"; transcript_id "T3";'],
            ['1', '.', 'CDS', '600', '649', '.', '+', '.',
             'gene_id "G2"; transcript_id "T3"; exon_number "1";'],
            ['1', '.', 'intron', '650', '749', '.', '+', '.',
             'gene_id "G2"; transcript_id "T3";'],
            ['1', '.', 'CDS', '750', '799', '.', '+', '.',
             'gene_id "G2"; transcript_id "T3"; exon_number "2";'],

        ])
        self.gtf = make_file_from_list(intervals_to_list(self.gtf_data))
        self.strange = get_temp_file_name()
        self.cross_tr = get_temp_file_name()
        self.out = get_temp_file_name()

    def test_explicit_whole_in(self):
        """
        Whole read is in single transcript and is crossing the exon-intron
        landmark (it is explicit). Provide three reads, with two different
        cross-links. One cross-link has two distinct randomers.
        """
        bam = make_bam_file({
            'chromosomes': [('1', 1000)],
            'segments': [
                # (qname, flag, refname, pos, mapq, cigar, tags)
                ('name2:rbc:CCCC', 0, 0, 140, 255, [(0, 50)], {'NH': 1}),
                ('name2:rbc:AAAA', 0, 0, 142, 255, [(0, 50)], {'NH': 1}),
                ('name2:rbc:CCCC', 0, 0, 142, 255, [(0, 50)], {'NH': 1}),
            ]
        })

        expected = [
            ['RNAmap', 'type', 'position', 'all', 'explicit'],
            ['UTR5-intron', '-10', '1', '1'],
            ['UTR5-intron', '-8', '2', '2'],
        ]

        rnamaps.run(bam, self.gtf, self.out, self.strange, self.cross_tr, mismatches=1)
        self.assertEqual(expected, make_list_from_file(self.out))

    def test_explicit_intergenic_left(self):
        """
        Read is half in intergenic region and half in transcript.
        """
        bam = make_bam_file({
            'chromosomes': [('1', 1000)],
            'segments': [
                # (qname, flag, refname, pos, mapq, cigar, tags)
                ('name2:rbc:CCCC', 0, 0, 80, 255, [(0, 50)], {'NH': 1}),
            ]
        })

        expected = [
            ['RNAmap', 'type', 'position', 'all', 'explicit'],
            ['intergenic-UTR5', '-20', '1', '1'],
        ]

        rnamaps.run(bam, self.gtf, self.out, self.strange, self.cross_tr, mismatches=1)
        self.assertEqual(expected, make_list_from_file(self.out))

    def test_explicit_intergenic_right(self):
        """
        Read is half in transcript region and half in intergenic.
        """
        bam = make_bam_file({
            'chromosomes': [('1', 1000)],
            'segments': [
                # (qname, flag, refname, pos, mapq, cigar, tags)
                ('name2:rbc:CCCC', 0, 0, 480, 255, [(0, 50)], {'NH': 1}),
            ]
        })

        expected = [
            ['RNAmap', 'type', 'position', 'all', 'explicit'],
            ['CDS-intergenic', '-20', '1', '1'],
        ]

        rnamaps.run(bam, self.gtf, self.out, self.strange, self.cross_tr, mismatches=1)
        self.assertEqual(expected, make_list_from_file(self.out))

    def test_cross_transcript_read(self):
        """
        Read is half in transcript region and half in intergenic.
        """
        bam = make_bam_file({
            'chromosomes': [('1', 1000)],
            'segments': [
                # (qname, flag, refname, pos, mapq, cigar, tags)
                ('name2:rbc:CCCC', 0, 0, 235, 255, [(0, 50)], {'NH': 1}),
            ]
        })

        expected = [
            ['chrom', 'strand', 'xlink', 'second-start', 'end-position', 'read_len'],
            ['1', '+', '234', '236', '284', '50'],
        ]

        rnamaps.run(bam, self.gtf, self.out, self.strange, self.cross_tr, mismatches=1)
        self.assertEqual(expected, make_list_from_file(self.cross_tr))

    def test_implicit_whole_in(self):
        """
        Whole read is in single transcript and in single segment. Also, this
        segment is the "middle" segment in transcript. Provide three reads, with
        two different cross-links. One cross-link has two distinct randomers.
        """
        bam = make_bam_file({
            'chromosomes': [('1', 1000)],
            'segments': [
                # (qname, flag, refname, pos, mapq, cigar, tags)
                ('name2:rbc:CCCC', 0, 0, 160, 255, [(0, 30)], {'NH': 1}),
                ('name2:rbc:CCCC', 0, 0, 163, 255, [(0, 30)], {'NH': 1}),
                ('name2:rbc:GGGG', 0, 0, 163, 255, [(0, 30)], {'NH': 1}),
            ]
        })

        expected = [
            ['RNAmap', 'type', 'position', 'all', 'explicit'],
            ['UTR5-intron', '10', '1', '0'],
            ['UTR5-intron', '13', '2', '0'],
        ]

        rnamaps.run(bam, self.gtf, self.out, self.strange, self.cross_tr, mismatches=1)
        self.assertEqual(expected, make_list_from_file(self.out))

    def test_implicit_exons(self):
        """
        Whole read is in single transcript and in single segment. Also, this
        segment is of EXON_TYPE in the "middle" segment in transcript. Only one read.
        """
        bam = make_bam_file({
            'chromosomes': [('1', 1000)],
            'segments': [
                # (qname, flag, refname, pos, mapq, cigar, tags)
                ('name2:rbc:CCCC', 0, 0, 205, 255, [(0, 20)], {'NH': 1}),
            ]
        })

        expected = [
            ['RNAmap', 'type', 'position', 'all', 'explicit'],
            ['CDS-UTR3', '-25', '0.25', '0'],
            ['CDS-intron', '-25', '0.25', '0'],
            ['UTR5-CDS', '5', '0.25', '0'],
            ['intron-CDS', '5', '0.25', '0'],
        ]

        rnamaps.run(bam, self.gtf, self.out, self.strange, self.cross_tr, mismatches=1,
                    implicit_handling='split')
        self.assertEqual(expected, make_list_from_file(self.out))

    def test_implicit_inter_tr(self):
        """
        Whole read is in single transcript, single segment. But the segment
        borders on intergenic (downstream).
        """
        bam = make_bam_file({
            'chromosomes': [('1', 1000)],
            'segments': [
                # (qname, flag, refname, pos, mapq, cigar, tags)
                ('name2:rbc:CCCC', 0, 0, 610, 255, [(0, 30)], {'NH': 1}),
            ]
        })

        expected = [
            ['RNAmap', 'type', 'position', 'all', 'explicit'],
            ['CDS-CDS', '-40', '0.3333', '0'],
            ['CDS-intron', '-40', '0.3333', '0'],
            ['intergenic-CDS', '10', '0.3333', '0'],
        ]

        rnamaps.run(bam, self.gtf, self.out, self.strange, self.cross_tr, mismatches=1,
                    implicit_handling='split')
        self.assertEqual(expected, make_list_from_file(self.out))

    def test_implicit_intergenic(self):
        """
        Whole read is in intergenic.
        """
        bam = make_bam_file({
            'chromosomes': [('1', 1000)],
            'segments': [
                # (qname, flag, refname, pos, mapq, cigar, tags)
                ('name2:rbc:CCCC', 0, 0, 530, 255, [(0, 30)], {'NH': 1}),
            ]
        })

        expected = [
            ['RNAmap', 'type', 'position', 'all', 'explicit'],
            ['CDS-intergenic', '30', '0.5', '0'],
            ['intergenic-CDS', '-70', '0.5', '0'],
        ]

        rnamaps.run(bam, self.gtf, self.out, self.strange, self.cross_tr, mismatches=1,
                    implicit_handling='split')
        self.assertEqual(expected, make_list_from_file(self.out))

    def test_negative_strand(self):
        """
        Whole read is in single transcript, single segment. But the segment
        borders on intergenic (downstream).
        """
        gtf_neg_data = [i[:6] + ['-'] + i[7:] for i in intervals_to_list(self.gtf_data)]
        gtf_neg = make_file_from_list(gtf_neg_data)
        bam = make_bam_file({
            'chromosomes': [('1', 1000)],
            'segments': [
                # (qname, flag, refname, pos, mapq, cigar, tags)
                ('name2:rbc:CCCC', 16, 0, 549, 255, [(0, 30)], {'NH': 1}),
            ]
        })

        expected = [
            ['RNAmap', 'type', 'position', 'all', 'explicit'],
            ['CDS-intergenic', '20', '0.5', '0'],
            ['intergenic-CDS', '-80', '0.5', '0'],
        ]

        rnamaps.run(bam, gtf_neg, self.out, self.strange, self.cross_tr, mismatches=1,
                    implicit_handling='split')
        self.assertEqual(expected, make_list_from_file(self.out))


class TestNormalisation(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", (ResourceWarning, ImportWarning))
        self.gtf_data = list_to_intervals([
            ['1', '.', 'intergenic', '1', '2', '.', '+', '.',
             'gene_id "."; transcript_id ".";'],
            # Gene #1:
            ['1', '.', 'gene', '3', '7', '.', '+', '.',
             'gene_id "G1";'],
            # Transcript #1
            ['1', '.', 'transcript', '3', '6', '.', '+', '.',
             'gene_id "G1"; transcript_id "T1";'],
            ['1', '.', 'CDS', '3', '3', '.', '+', '.',
             'gene_id "G1"; transcript_id "T1"; exon_number "2";'],
            ['1', '.', 'intron', '4', '6', '.', '+', '.',
             'gene_id "G1"; transcript_id "T1";'],
            ['1', '.', 'UTR3', '5', '6', '.', '+', '.',
             'gene_id "G1"; transcript_id "T1"; exon_number "3";'],

            # Transcript #2
            ['1', '.', 'transcript', '4', '7', '.', '+', '.',
             'gene_id "G1"; transcript_id "T2";'],
            ['1', '.', 'ncRNA', '4', '5', '.', '+', '.',
             'gene_id "G1"; transcript_id "T2"; exon_number "1";'],
            ['1', '.', 'intron', '6', '6', '.', '+', '.',
             'gene_id "G1"; transcript_id "T2";'],
            ['1', '.', 'ncRNA', '7', '7', '.', '+', '.',
             'gene_id "G1"; transcript_id "T2"; exon_number "2";'],

            # intergenic
            ['1', '.', 'intergenic', '8', '9', '.', '+', '.',
             'gene_id "."; transcript_id ".";'],
        ])
        self.gtf = make_file_from_list(intervals_to_list(self.gtf_data))

    def test_normalisation(self):
        norm_file = get_temp_file_name(extension='txt')
        rnamaps.make_normalization(self.gtf, norm_file)

        expected = [
            ['RNAmap_type', 'distance', 'segments'],
            ['CDS-UTR3', '-1', '1'],
            ['CDS-UTR3', '0', '1'],
            ['CDS-UTR3', '1', '1'],
            ['CDS-intron', '-1', '1'],
            ['CDS-intron', '0', '1'],
            ['CDS-intron', '1', '1'],
            ['CDS-intron', '2', '1'],
            ['integrenic-CDS', '-2', '1'],
            ['integrenic-CDS', '-1', '1'],
            ['integrenic-CDS', '0', '1'],
            ['intron-UTR3', '-3', '1'],
            ['intron-UTR3', '-2', '1'],
            ['intron-UTR3', '-1', '1'],
            ['intron-UTR3', '0', '1'],
            ['intron-UTR3', '1', '1'],
            ['intron-ncRNA', '-1', '1'],
            ['intron-ncRNA', '0', '1'],
            ['ncRNA-integrenic', '-1', '1'],
            ['ncRNA-integrenic', '0', '1'],
            ['ncRNA-integrenic', '1', '1'],
            ['ncRNA-intron', '-2', '1'],
            ['ncRNA-intron', '-1', '1'],
            ['ncRNA-intron', '0', '1'],
            ['ncRNA-ncRNA', '-2', '1'],
            ['ncRNA-ncRNA', '-1', '1'],
            ['ncRNA-ncRNA', '0', '1'],
        ]

        self.assertEqual(expected, make_list_from_file(norm_file))

    def test_plot(self):
        image_file = get_temp_file_name(extension='png')
        norm_file = get_temp_file_name(extension='txt')
        rnamaps.make_normalization(self.gtf, norm_file)
        rnamaps.plot_rna_map(norm_file, 'CDS-intron', normalization=norm_file, outfile=image_file)
        self.assertTrue(os.path.isfile(image_file))


if __name__ == '__main__':
    unittest.main()
