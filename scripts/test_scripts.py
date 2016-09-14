"""
This script tests if CLI interface works as expexcted.
"""

import os
import unittest
import subprocess

from iCount.tests.utils import make_bam_file, get_temp_file_name, \
    make_file_from_list, get_temp_dir


class TestCLI(unittest.TestCase):

    def setUp(self):
        # Temporary file names to use for output:
        self.tmp1 = get_temp_file_name()
        self.tmp2 = get_temp_file_name()
        self.dir = get_temp_dir()

        # Make a sample bam file
        self.bam_data = {
            'chromosomes': [('chr1', 3000), ('chr2', 2000)],
            'segments': [
                ('name1', 4, 0, 100, 20, [(0, 100)], {'NH': 1}),
                ('name2', 0, 0, 100, 20, [(0, 201)], {'NH': 7}),
                ('name3:rbc:CCCC:', 16, 0, 100, 20, [(0, 100)], {'NH': 1}),
                ('name4:ABC', 0, 1, 300, 20, [(0, 200)], {'NH': 11}),
                ('name4:ACG', 0, 1, 300, 20, [(0, 200)], {'NH': 11}),
                ('name5', 0, 1, 300, 3, [(0, 200)], {'NH': 13})]}
        self.bam = make_bam_file(self.bam_data)

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

        self.chrom_len = make_file_from_list(bedtool=False, data=[
            ['1', '2000'],
            ['2', '1000'],
        ])

    def test_annotate(self):
        command_basic = ['iCount', 'annotate', self.annotation,
                         self.cross_links, self.tmp1]
        command_full = ['iCount', 'annotate', self.annotation,
                        self.cross_links, self.tmp1,
                        '--subtype', 'biotype',
                        '--excluded_types', 'ncRNA']

        self.assertEqual(subprocess.call(command_basic), 0)
        self.assertEqual(subprocess.call(command_full), 0)

    def test_clusters(self):
        command_basic = ['iCount', 'clusters', self.cross_links,
                         self.tmp1]
        command_full = ['iCount', 'clusters', self.cross_links,
                        self.tmp1, '--dist', '20']

        self.assertEqual(subprocess.call(command_basic), 0)
        self.assertEqual(subprocess.call(command_full), 0)

    def test_group(self):
        command_basic = ['iCount', 'group', self.tmp1,
                         # Simulate list of two files:
                         self.cross_links + ', ', self.cross_links]

        self.assertEqual(subprocess.call(command_basic), 0)

    def test_peaks(self):
        command_basic = ['iCount', 'peaks', self.annotation,
                         self.cross_links, self.tmp1]
        command_full = [
            'iCount', 'peaks', self.annotation,
            self.cross_links, self.tmp1,
            '--fout_scores', self.tmp2,
            '--hw', '3',
            '--fdr', '0.05',
            '--perms', '10',
            '--rnd_seed', '42',
            '--features', 'gene',
            '--report_progress'
        ]

        self.assertEqual(subprocess.call(command_basic), 0)
        self.assertEqual(subprocess.call(command_full), 0)

    def test_summary(self):
        command_basic = ['iCount', 'summary', self.annotation,
                         self.cross_links, self.tmp1, self.chrom_len]
        command_full = ['iCount', 'summary', self.annotation,
                        self.cross_links, self.tmp1, self.chrom_len,
                        '--ndigits', '8',
                        '--subtype', 'biotype',
                        '--excluded_types', 'ncRNA,', 'intron']

        self.assertEqual(subprocess.call(command_basic), 0)
        self.assertEqual(subprocess.call(command_full), 0)

    # #################################
    # #################################

    def test_releases(self):
        # TODO: rename to get_releases?
        command_basic = ['iCount', 'get_release_list']
        self.assertEqual(subprocess.call(command_basic), 0)

    def test_species(self):
        # TODO: rename to get_species?
        command_basic = ['iCount', 'get_species_list', '84']
        self.assertEqual(subprocess.call(command_basic), 0)

    def test_get_annotation(self):
        # TODO: rename to get_annotation?
        command_basic = ['iCount', 'download_annotation', '84', 'homo_sapiens']
        command_full = ['iCount', 'download_annotation', '84', 'homo_sapiens',
                        '--target_dir', self.dir,
                        '--target_fname', self.tmp1,
                        ]

        self.assertEqual(subprocess.call(command_basic), 0)
        self.assertEqual(subprocess.call(command_full), 0)

    def test_get_sequence(self):
        # TODO: rename to get_genome?

        # Download just MT chromosome, or this can last too long...
        command_full = ['iCount', 'download_sequence', '84', 'homo_sapiens',
                        '--target_dir', self.dir,
                        '--target_fname', self.tmp1,
                        '--tempdir', self.dir,
                        '--chromosomes', 'MT', 'Y',
                        ]

        self.assertEqual(subprocess.call(command_full), 0)

    # #################################
    # #################################

    # TODO: make this test!
    @unittest.skip
    def test_map(self):
        # self.reads
        # self.dir = tempdir...
        command_basic = ['iCount', 'map', '-h']

        self.assertEqual(subprocess.call(command_basic), 0)
        self.assertEqual(subprocess.call(command_full), 0)

    # TODO: make this test!
    @unittest.skip
    def test_mapindex(self):
        command_basic = ['iCount', 'mapindex', '-h']

        self.assertEqual(subprocess.call(command_basic), 0)
        self.assertEqual(subprocess.call(command_full), 0)

    def test_xlsites(self):
        command_basic = ['iCount', 'xlsites', self.bam, self.tmp1, self.tmp2]
        command_full = ['iCount', 'xlsites', self.bam, self.tmp1, self.tmp2,
                        '--group_by', 'start',
                        '--quant', 'cDNA',
                        '--randomer_mismatches', '2',
                        '--mapq_th', '0',
                        '--multimax', '50'
                        ]

        self.assertEqual(subprocess.call(command_basic), 0)
        self.assertEqual(subprocess.call(command_full), 0)


if __name__ == '__main__':
    unittest.main()
