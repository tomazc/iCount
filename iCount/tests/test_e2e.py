# pylint: disable=missing-docstring, protected-access
import os
import unittest
import warnings

from numpy import random

import iCount
from iCount.tests.utils import get_temp_file_name, get_temp_dir, make_file_from_list, \
    make_fasta_file, make_sequence, make_quality_scores, attrs


class TestEndToEnd(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", (ResourceWarning, ImportWarning))

        # Make fai and fasta file:
        self.fai_data = {'1': 1000}
        self.fai = make_file_from_list(list(self.fai_data.items()), bedtool=False, extension='fai')

        random.seed(0)  # pylint:disable=no-member
        random_seeds = random.randint(10**5, size=len(self.fai))  # pylint:disable=no-member
        self.seqs = {name: make_sequence(size, rnd_seed=rseed) for (name, size), rseed in
                     zip(self.fai_data.items(), random_seeds)}
        self.fasta = make_fasta_file(sequences=self.seqs.values(), headers=self.seqs.keys())

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
            ['1', '.', 'gene', '600', '899', '.', '+', '.', attrs(gid='G2')],
            # Transcript #3
            ['1', '.', 'transcript', '600', '899', '.', '+', '.', attrs(gid='G2', tid='T3')],
            ['1', '.', 'exon', '600', '659', '.', '+', '.', attrs(gid='G2', tid='T3', exn=1)],
            ['1', '.', 'CDS', '600', '659', '.', '+', '.', attrs(gid='G2', tid='T3', exn=1)],
            ['1', '.', 'exon', '840', '899', '.', '+', '.', attrs(gid='G2', tid='T3', exn=2)],
            ['1', '.', 'CDS', '840', '899', '.', '+', '.', attrs(gid='G2', tid='T3', exn=2)],

            # Gene #3, negative strand:
            ['1', '.', 'gene', '800', '899', '.', '-', '.', 'gene_id "G3";'],
            # Transcript #4
            ['1', '.', 'transcript', '800', '899', '.', '-', '.', attrs(gid='G3', tid='T4')],
            ['1', '.', 'exon', '800', '899', '.', '-', '.', attrs(gid='G3', tid='T4', exn=1)],
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
        random_seeds = random.randint(10**5, size=len(self.read_data))  # pylint:disable=no-member
        self.reads = get_temp_file_name(extension='fastq')
        with open(self.reads, 'wt') as ofile:
            for (name, chrom, length, xlink, strand), rseed in zip(self.read_data, random_seeds):
                seq = self.seqs[chrom][xlink: xlink + length]
                if strand == '-':
                    seq = self.seqs[chrom][xlink - length: xlink]
                    compl_bases = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
                    seq = ''.join([compl_bases[nuc] for nuc in seq[::-1]])

                ofile.write('@' + name + '\n')
                ofile.write(seq + '\n')
                ofile.write('+' + '\n')
                ofile.write(
                    make_quality_scores(len(seq), min_chr=65, max_chr=73, rnd_seed=rseed) + '\n')

            # Make a "strange" read:
            # pylint: disable=undefined-loop-variable
            seq = self.seqs[chrom][250: 295] + self.seqs[chrom][310: 380]
            ofile.write('@' + 'name_strange:rbc:GGGG' + '\n')
            ofile.write(seq + '\n')
            ofile.write('+' + '\n')
            ofile.write(
                make_quality_scores(len(seq), min_chr=70, max_chr=73, rnd_seed=rseed) + '\n')

    def test_e2e(self):
        """
        From raw reads and ENSEMBL annotation to rnamaps.
        """

        # Make segmentation & regions file
        seg = get_temp_file_name(extension='gtf')
        out_dir = get_temp_dir()
        iCount.genomes.segment.get_segments(self.gtf, seg, self.fai)
        iCount.genomes.region.make_regions(seg, out_dir)
        regions = os.path.join(out_dir, iCount.genomes.region.REGIONS_FILE)

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

        sites_single = get_temp_file_name(extension='bed.gz')
        sites_multi = get_temp_file_name(extension='bed.gz')
        skipped = get_temp_file_name(extension='bam')
        iCount.mapping.xlsites.run(bam, sites_single, sites_multi, skipped)

        iCount.analysis.rnamaps.run(sites_single, regions)


if __name__ == '__main__':
    unittest.main()
