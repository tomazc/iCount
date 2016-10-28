"""
This script tests if CLI interface works as expexcted.
"""

import os
import unittest
import subprocess
import warnings

from iCount.tests.utils import make_bam_file, get_temp_file_name, \
    make_file_from_list, get_temp_dir, make_fastq_file, make_fasta_file


class TestCLI(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", ResourceWarning)

        # Temporary file names to use for output:
        self.tmp1 = get_temp_file_name()
        self.tmp2 = get_temp_file_name()
        self.dir = get_temp_dir()
        self.dir2 = get_temp_dir()

        self.cross_links = make_file_from_list([
            ['1', '16', '17', '.', '5', '+'],
            ['1', '14', '15', '.', '5', '+'],
            ['1', '15', '16', '.', '5', '+'],
        ])

        self.annotation = make_file_from_list([
            ['1', '.', 'CDS', '10', '20', '.', '+', '.', 'biotype "A";'],
            ['1', '.', 'ncRNA', '10', '20', '.', '+', '.', 'biotype "A";'],
            ['1', '.', 'CDS', '10', '20', '.', '+', '.', 'biotype "A";'],
            ['1', '.', 'CDS', '10', '20', '.', '+', '.', 'biotype "B";'],
            ['1', '.', 'CDS', '10', '20', '.', '-', '.', 'biotype "C";'],
            ['1', '.', 'CDS', '12', '18', '.', '+', '.', 'biotype "A";'],
            ['1', '.', 'CDS', '30', '40', '.', '+', '.', 'biotype "D";'],
        ])

        self.gtf = make_file_from_list([
            ['1', '.', 'gene', '10', '20', '.', '+', '.',
             'gene_id "A";'],
            ['1', '.', 'transcript', '10', '20', '.', '+', '.',
             'gene_id "A"; transcript_id "AA";'],
            ['1', '.', 'exon', '10', '20', '.', '+', '.',
             'gene_id "A"; transcript_id "AA"; exon_number "1";'],
        ])

    def test_releases(self):
        command_basic = ['iCount', 'releases',
                         '-S', '40',  # Supress lower than ERROR messages.
                         ]
        self.assertEqual(subprocess.call(command_basic), 0)

    def test_species(self):
        command_basic = ['iCount', 'species',
                         '-S', '40',  # Supress lower than ERROR messages.
                         ]
        self.assertEqual(subprocess.call(command_basic), 0)

    def test_annotation(self):
        # Execute only full command (with --target_dir), to avoid downloading to cwd.
        command_full = ['iCount', 'annotation', 'homo_sapiens',
                        '--release', '84',
                        '--out_dir', self.dir,
                        '--annotation', get_temp_file_name(extension='gtf.gz'),
                        '-S', '40',  # Supress lower than ERROR messages.
                        ]

        self.assertEqual(subprocess.call(command_full), 0)

    def test_genome(self):
        # Download just MT and Y chromosome, or test can last too long...
        command_full = ['iCount', 'genome', 'homo_sapiens',
                        '--release', '84',
                        '--out_dir', self.dir,
                        '--genome', get_temp_file_name(extension='fa.gz'),
                        '--chromosomes', 'MT', 'Y',
                        '-S', '40',  # Supress lower than ERROR messages.
                        ]

        self.assertEqual(subprocess.call(command_full), 0)

    def test_genes(self):
        fai = make_file_from_list([
                ['1', '2000'],
                ['MT', '500'],
        ], bedtool=False)

        command_basic = ['iCount', 'genes', self.gtf, self.tmp1,
                         '-S', '40',  # Supress lower than ERROR messages.
                         ]
        command_full = ['iCount', 'genes', self.gtf, self.tmp1,
                        '--fai', fai,
                        '--feature', 'gene',
                        '--attribute', 'gene_id',
                        '-S', '40',  # Supress lower than ERROR messages.
                        ]

        self.assertEqual(subprocess.call(command_basic), 0)
        self.assertEqual(subprocess.call(command_full), 0)

    def test_segment(self):
        fai = make_file_from_list([
                ['1', '2000'],
                ['MT', '500'],
        ], bedtool=False)

        command_basic = ['iCount', 'segment', self.gtf, self.tmp1, fai,
                         '-S', '40',  # Supress lower than ERROR messages.
                         ]
        command_full = ['iCount', 'segment', self.gtf, self.tmp1, fai,
                        '--report_progress',
                        '-S', '40',  # Supress lower than ERROR messages.
                        ]

        self.assertEqual(subprocess.call(command_basic), 0)
        self.assertEqual(subprocess.call(command_full), 0)

    # #################################
    # #################################

    def test_demultiplex(self):
        barcodes = [
            'NNNTTGTNN',
            'NNNGGTTNN',
            'NNNGGCGNN',
        ]
        b1, b2, b3 = barcodes
        adapter = 'CCCCCCCCC'
        fastq = make_fastq_file(barcodes=barcodes, adapter=adapter)

        command_basic = ['iCount', 'demultiplex',
                         fastq,
                         adapter,
                         b1, b2, b3,
                         '--mismatches', '2',
                         '--out_dir', self.dir,  # put files in tmpdir, to not pollute cwd.
                         '-S', '40',  # Supress lower than ERROR messages.
                         ]
        command_full = ['iCount', 'demultiplex',
                        fastq,
                        adapter,
                        b1, b2, b3,
                        '--mismatches', '2',
                        '--minimum_length', '15',
                        '--prefix', 'demux',
                        '--out_dir', self.dir,
                        '-S', '40',  # Supress lower than ERROR messages.
                        ]

        self.assertEqual(subprocess.call(command_basic), 0)
        self.assertEqual(subprocess.call(command_full), 0)

    def test_cutadapt(self):
        adapter = 'CCCCCCCCC'
        fastq = make_fastq_file(adapter=adapter, out_file=get_temp_file_name(extension='fastq'))

        command_basic = ['iCount', 'cutadapt',
                         fastq,
                         self.tmp1,
                         adapter,
                         '-S', '40',  # Supress lower than ERROR messages.
                         ]
        command_full = ['iCount', 'cutadapt',
                        fastq,
                        self.tmp1,
                        adapter,
                        '--qual_base', '64',
                        '--qual_trim', '20',
                        '--minimum_length', '15',
                        '-S', '40',  # Supress lower than ERROR messages.
                        ]

        self.assertEqual(subprocess.call(command_basic), 0)
        self.assertEqual(subprocess.call(command_full), 0)

    def test_mapindex_and_map_basic(self):
        genome = make_fasta_file(num_sequences=2, seq_len=1000)
        fastq = make_fastq_file(genome=genome)

        command_basic = ['iCount', 'mapindex', genome, self.dir,
                         '-S', '40',  # Supress lower than ERROR messages.
                         ]
        self.assertEqual(subprocess.call(command_basic), 0)

        command_basic = ['iCount', 'map',
                         fastq,
                         self.dir,
                         self.dir2,
                         '-S', '40',  # Supress lower than ERROR messages.
                         ]

        self.assertEqual(subprocess.call(command_basic), 0)

    def test_mapindex_and_map_full(self):
        genome = make_fasta_file(num_sequences=2, seq_len=1000)
        fastq = make_fastq_file(genome=genome)

        command_full = ['iCount', 'mapindex', genome, self.dir,
                        '--annotation', self.gtf,
                        '--overhang', '100',
                        '--overhang_min', '8',
                        '--threads', '1',
                        '-S', '40',  # Supress lower than ERROR messages.
                        ]
        self.assertEqual(subprocess.call(command_full), 0)

        command_full = ['iCount', 'map',
                        fastq,
                        self.dir,
                        self.dir2,
                        '--annotation', self.gtf,
                        '--multimax', '50',
                        '--mismatches', '2',
                        '--threads', '1',
                        '-S', '40',  # Supress lower than ERROR messages.
                        ]
        self.assertEqual(subprocess.call(command_full), 0)

    def test_xlsites(self):
        # Make a sample bam file
        unique = get_temp_file_name(extension='.bed')
        multi = get_temp_file_name(extension='.bed')
        strange = get_temp_file_name(extension='.bam')
        bam = make_bam_file({
            'chromosomes': [
                ('chr1', 3000),
                ('chr2', 2000),
            ],
            'segments': [
                ('name1', 4, 0, 100, 20, [(0, 100)], {'NH': 1}),
                ('name2', 0, 0, 100, 20, [(0, 201)], {'NH': 7}),
                ('name3:rbc:CCCC:', 16, 0, 100, 20, [(0, 100)], {'NH': 1}),
                ('name4:ABC', 0, 1, 300, 20, [(0, 200)], {'NH': 11}),
                ('name4:ACG', 0, 1, 300, 20, [(0, 200)], {'NH': 11}),
                ('name5', 0, 1, 300, 3, [(0, 200)], {'NH': 13})],
            })

        command_basic = ['iCount', 'xlsites', bam, unique, multi, strange,
                         '-S', '40',  # Supress lower than ERROR messages.
                         ]
        command_full = ['iCount', 'xlsites', bam, unique, multi, strange,
                        '--group_by', 'start',
                        '--quant', 'cDNA',
                        '--mismatches', '2',
                        '--mapq_th', '0',
                        '--multimax', '50',
                        '-S', '40',  # Supress lower than ERROR messages.
                        ]

        self.assertEqual(subprocess.call(command_basic), 0)
        self.assertEqual(subprocess.call(command_full), 0)

    # #################################
    # #################################

    def test_annotate(self):
        command_basic = ['iCount', 'annotate', self.annotation,
                         self.cross_links, self.tmp1,
                         '-S', '40',  # Supress lower than ERROR messages.
                         ]
        command_full = ['iCount', 'annotate', self.annotation,
                        self.cross_links, self.tmp1,
                        '--subtype', 'biotype',
                        '--excluded_types', 'ncRNA',
                        '-S', '40',  # Supress lower than ERROR messages.
                        ]

        self.assertEqual(subprocess.call(command_basic), 0)
        self.assertEqual(subprocess.call(command_full), 0)

    def test_clusters(self):
        command_basic = ['iCount', 'clusters', self.cross_links,
                         self.tmp1,
                         '-S', '40',  # Supress lower than ERROR messages.
                         ]
        command_full = ['iCount', 'clusters', self.cross_links,
                        self.tmp1, '--dist', '20',
                        '-S', '40',  # Supress lower than ERROR messages.
                        ]

        self.assertEqual(subprocess.call(command_basic), 0)
        self.assertEqual(subprocess.call(command_full), 0)

    def test_group(self):
        command_basic = ['iCount', 'group', self.tmp1,
                         # Simulate list of two files:
                         self.cross_links + ', ', self.cross_links,
                         '-S', '40',  # Supress lower than ERROR messages.
                         ]

        self.assertEqual(subprocess.call(command_basic), 0)

    def test_peaks(self):
        command_basic = ['iCount', 'peaks', self.annotation,
                         self.cross_links, get_temp_file_name(extension='.bed.gz'),
                         '-S', '40',  # Supress lower than ERROR messages.
                         ]
        command_full = [
            'iCount', 'peaks', self.annotation,
            self.cross_links, get_temp_file_name(extension='.bed.gz'),
            '--scores', get_temp_file_name(extension='.tsv.gz'),
            '--hw', '3',
            '--fdr', '0.05',
            '--perms', '10',
            '--rnd_seed', '42',
            '--features', 'gene',
            '--report_progress',
            '-S', '40',  # Supress lower than ERROR messages.
        ]

        self.assertEqual(subprocess.call(command_basic), 0)
        self.assertEqual(subprocess.call(command_full), 0)

    def test_summary(self):
        chrom_len = make_file_from_list(bedtool=False, data=[
            ['1', '2000'],
            ['2', '1000'],
        ])

        command_basic = ['iCount', 'summary', self.annotation,
                         self.cross_links, self.tmp1, chrom_len,
                         '-S', '40',  # Supress lower than ERROR messages.
                         ]
        command_full = ['iCount', 'summary', self.annotation,
                        self.cross_links, self.tmp1, chrom_len,
                        '--digits', '8',
                        '--subtype', 'biotype',
                        '--excluded_types', 'ncRNA,', 'intron',
                        '-S', '40',  # Supress lower than ERROR messages.
                        ]

        self.assertEqual(subprocess.call(command_basic), 0)
        self.assertEqual(subprocess.call(command_full), 0)

if __name__ == '__main__':
    unittest.main()
