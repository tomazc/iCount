# pylint: disable=missing-docstring, protected-access
import os
import unittest
import warnings

import pysam

import iCount
from iCount.tests.utils import get_temp_file_name, get_temp_dir, make_file_from_list, \
    make_list_from_file, make_fasta_file, make_sequence, make_quality_scores, attrs


class TestEndToEnd(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", (ResourceWarning, ImportWarning))

        # Make fai and fasta file:
        self.fai_data = {'1': 1000}
        self.fai = make_file_from_list(list(self.fai_data.items()), bedtool=False, extension='fai')

        self.seqs = {name: make_sequence(size) for name, size in self.fai_data.items()}
        self.fasta = make_fasta_file(sequences=self.seqs.values(),
                                     headers=self.seqs.keys())

        # Make annotation:
        self.gtf = make_file_from_list(extension='gtf', data=[
            # Gene #1, positive strand:
            ['1', '.', 'gene', '100', '499', '.', '+', '.', attrs(gid='G1')],
            # Transcript #1
            ['1', '.', 'transcript', '100', '249', '.', '+', '.', attrs(gid='G1', tid='T1')],
            ['1', '.', 'exon', '100', '149', '.', '+', '.', attrs(gid='G1', tid='T1', exn=1)],
            ['1', '.', 'exon', '200', '229', '.', '+', '.', attrs(gid='G1', tid='T1', exn=2)],
            ['1', '.', 'CDS', '200', '229', '.', '+', '.', attrs(gid='G1', tid='T1', exn=2)],
            ['1', '.', 'exon', '240', '249', '.', '+', '.', attrs(gid='G1', tid='T1', exn=3)],
            # Transcript #2
            ['1', '.', 'transcript', '240', '499', '.', '+', '.', attrs(gid='G1', tid='T2')],
            ['1', '.', 'exon', '240', '299', '.', '+', '.', attrs(gid='G1', tid='T2', exn=1)],
            ['1', '.', 'CDS', '240', '299', '.', '+', '.', attrs(gid='G1', tid='T2', exn=1)],
            ['1', '.', 'exon', '400', '499', '.', '+', '.', attrs(gid='G1', tid='T2', exn=2)],
            ['1', '.', 'CDS', '400', '499', '.', '+', '.', attrs(gid='G1', tid='T2', exn=2)],

            # Gene #2, positive strand:
            ['1', '.', 'gene', '600', '799', '.', '+', '.', attrs(gid='G2')],
            # Transcript #3
            ['1', '.', 'transcript', '600', '799', '.', '+', '.', attrs(gid='G2', tid='T3')],
            ['1', '.', 'exon', '600', '649', '.', '+', '.', attrs(gid='G2', tid='T3', exn=1)],
            ['1', '.', 'CDS', '600', '649', '.', '+', '.', attrs(gid='G2', tid='T3', exn=1)],
            ['1', '.', 'exon', '750', '799', '.', '+', '.', attrs(gid='G2', tid='T3', exn=2)],
            ['1', '.', 'CDS', '750', '799', '.', '+', '.', attrs(gid='G2', tid='T3', exn=2)],

            # Gene #3, negative strand:
            ['1', '.', 'gene', '800', '899', '.', '-', '.', 'gene_id "G3";'],
            # Transcript #4
            ['1', '.', 'transcript', '800', '899', '.', '-', '.', attrs(gid='G3', tid='T4')],
            ['1', '.', 'exon', '800', '899', '.', '-', '.', attrs(gid='G3', tid='T4', exn=1)],
            ['1', '.', 'CDS', '800', '899', '.', '-', '.', attrs(gid='G3', tid='T4', exn=1)],
        ])

        # Define positions of cross-links:
        self.read_data = [
            # [name, chromosome, sequence_len, read_start, strand]
            ('name01:rbc:CCCC', '1', 50, 80, '+'),
            ('name02:rbc:CCCC', '1', 50, 140, '+'),
            ('name03:rbc:AAAA', '1', 50, 142, '+'),
            ('name04:rbc:CCCC', '1', 50, 142, '+'),
            ('name05:rbc:CCCC', '1', 30, 160, '+'),
            ('name06:rbc:CCCC', '1', 30, 163, '+'),
            ('name07:rbc:GGGG', '1', 30, 163, '+'),
            ('name08:rbc:CCCC', '1', 20, 205, '+'),
            ('name09:rbc:CCCC', '1', 50, 235, '+'),
            ('name10:rbc:CCCC', '1', 50, 480, '+'),
            ('name11:rbc:CCCC', '1', 30, 530, '+'),
            ('name12:rbc:CCCC', '1', 30, 610, '+'),
            ('name13:rbc:CCCC', '1', 30, 549, '-'),
        ]
        # Make reads from genome, pointing to specific pre-detemrined cross-link positions:
        self.reads = get_temp_file_name(extension='fastq')
        with open(self.reads, 'wt') as ofile:
            for name, chrom, length, xlink, strand in self.read_data:
                seq = self.seqs[chrom][xlink: xlink + length]
                if strand == '-':
                    seq = self.seqs[chrom][xlink - length: xlink]
                    compl_bases = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
                    seq = ''.join([compl_bases[nuc] for nuc in seq[::-1]])

                ofile.write('@' + name + '\n')
                ofile.write(seq + '\n')
                ofile.write('+' + '\n')
                ofile.write(make_quality_scores(len(seq), min_chr=65, max_chr=73) + '\n')

            # Make a "strange" read:
            # pylint: disable=undefined-loop-variable
            seq = self.seqs[chrom][250: 295] + self.seqs[chrom][310: 380]
            ofile.write('@' + 'name_strange:rbc:GGGG' + '\n')
            ofile.write(seq + '\n')
            ofile.write('+' + '\n')
            ofile.write(make_quality_scores(len(seq), min_chr=70, max_chr=73) + '\n')

    def test_e2e(self):
        """
        From raw reads and ENSEMBL annotation to rnamaps.
        """

        # Make segmentation
        seg = get_temp_file_name(extension='gtf')
        iCount.genomes.segment.get_regions(self.gtf, seg, self.fai)

        # Build STAR index:
        genome_index = get_temp_dir()
        rcode = iCount.externals.star.build_index(self.fasta, genome_index, annotation=self.gtf)
        self.assertEqual(rcode, 0)
        # Map reads:
        map_dir = get_temp_dir()
        rcode = iCount.externals.star.map_reads(
            self.reads, genome_index, out_dir=map_dir, annotation=self.gtf)
        self.assertEqual(rcode, 0)

        # Get bam with mapped reads:
        bam = [fname for fname in os.listdir(map_dir) if fname.startswith('Aligned')][0]
        bam = os.path.join(map_dir, bam)

        # Make all sorts of analysis and save it:
        normal_out = get_temp_file_name(extension='tsv')
        strange_out = get_temp_file_name(extension='bam')
        cross_tr_out = get_temp_file_name(extension='tsv')
        iCount.analysis.rnamaps.run(bam, seg, normal_out, strange_out, cross_tr_out,
                                    implicit_handling='split')
        # Normal output:
        expected_out = [
            ['RNAmap', 'type', 'position', 'all', 'explicit'],
            ['CDS-CDS', '-40', '0.3333', '0'],
            ['CDS-UTR3', '-25', '0.25', '0'],
            ['CDS-intergenic', '-20', '1', '1'],
            ['CDS-intergenic', '30', '0.5', '0'],
            ['CDS-intergenic', '250', '1', '0'],
            ['CDS-intron', '-40', '0.3333', '0'],
            ['CDS-intron', '-25', '0.25', '0'],
            ['UTR5-CDS', '5', '0.25', '0'],
            ['UTR5-intron', '-10', '1', '1'],
            ['UTR5-intron', '-8', '2', '2'],
            ['UTR5-intron', '10', '0.5', '0'],
            ['UTR5-intron', '13', '1', '0'],
            ['intergenic-CDS', '-70', '0.5', '0'],
            ['intergenic-CDS', '10', '0.3333', '0'],
            ['intergenic-UTR5', '-20', '1', '1'],
            ['intron-CDS', '-40', '0.5', '0'],
            ['intron-CDS', '-37', '1', '0'],
            ['intron-CDS', '5', '0.25', '0'],
        ]
        self.assertEqual(expected_out, make_list_from_file(normal_out))

        # Cross transcript:
        expected_cross_transcript = [
            ['chrom', 'strand', 'xlink', 'second-start', 'end-position', 'read_len'],
            ['1', '+', '234', '236', '284', '50'],
        ]
        self.assertEqual(expected_cross_transcript, make_list_from_file(cross_tr_out))

        # Strange:
        strange_reads = list(pysam.AlignmentFile(strange_out, 'rb'))  # pylint: disable=no-member
        self.assertEqual(len(strange_reads), 1)
        strange_read = strange_reads[0]
        self.assertEqual(strange_read.query_name, 'name_strange:rbc:GGGG')
        self.assertEqual(strange_read.reference_start, 250)
        self.assertEqual(strange_read.cigar, [(0, 45), (2, 15), (0, 70)])


if __name__ == '__main__':
    unittest.main()
