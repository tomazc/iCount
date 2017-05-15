"""
Test iCount.genomes.segments.get_regions consistency.

This script tests if results of iCount.genomes.segments.get_regions function
are consistent with the reference results. If reference results are not
present, they are created.

If gtf/fasta/genome/get_regions.gtf files for given release and species are
not found in `self.files_dir` they are downloaded/calculated. This may
significantly extend the execution time of this script.

By modifying the variables `self.releases_list` and `self.species_list` user can
determine which gtf files to put under test.
"""
# pylint: disable=missing-docstring, protected-access

import unittest
import os
from os.path import join, exists, isfile
from uuid import uuid4
from collections import Counter

import pybedtools

from iCount.genomes import ensembl, segment
from iCount.tests.functional import FUNCTIONAL_TEST_FOLDER


class TestSegmentation(unittest.TestCase):

    def setUp(self):
        self.ref_results_file = __file__[:-3] + '.out.txt'

        # Folder where all gtf files are stored:
        self.files_dir = join(FUNCTIONAL_TEST_FOLDER, 'gtf_files/')
        if not exists(self.files_dir):
            os.makedirs(self.files_dir)

        self.releases_list = ensembl.releases()
        self.species_list = ['homo_sapiens', 'mus_musculus']

        self.order = ['gene', 'transcript', 'CDS', 'UTR3', 'UTR5', 'intron',
                      'ncRNA', 'intergenic']
        header = '\t'.join(map(str, ['species', 'release'] + [o for o in self.order]))

        self.current_results = [header]
        self.reference_results = []
        self.success = False

    def tearDown(self):
        if not self.success and self.current_results:
            if isfile(self.ref_results_file):
                self.ref_results_file = self.ref_results_file + str(uuid4())
            with open(self.ref_results_file, 'w') as ref_file:
                for result in self.current_results:
                    ref_file.write(result + '\n')

    def cacluate_gtf_regions(self, species, release):
        """
        Calculate gtf_regions for given release and species.

        Also, download/compute all the necessary files.
        """
        gtf_regions = join(self.files_dir, '{}.{}.get_regions.gtf.gz'.format(species, release))

        if not isfile(gtf_regions):
            # If gtf_regions file is not preswent, make it from gtf_base and chrom_length file

            gtf_base = join(self.files_dir, '{}.{}.gtf.gz'.format(species, release))
            if not isfile(gtf_base):
                ensembl.annotation(species, release=release, out_dir=self.files_dir)

            chrom_length = join(
                self.files_dir, '{}.{}.fa.gz.chrom_length.txt'.format(species, release))
            if not isfile(chrom_length):
                fasta = join(self.files_dir, '{}.{}.fa.gz'.format(species, release))
                if not isfile(fasta):
                    fasta = ensembl.genome(release, species, out_dir=self.files_dir)
                # Now we have all files to produce chrom_length
                chrom_length = ensembl.chrom_length(fasta)

            # Now we have gtf_base and chrom_length: produce gtf_regions:
            gtf_regions = segment.get_regions(gtf_base, gtf_regions, chrom_length)

        return gtf_regions

    def test_reference_stats(self):
        for species in self.species_list:
            for release in self.releases_list:
                gtf_regions = self.cacluate_gtf_regions(species, release)
                result = Counter(i[2] for i in pybedtools.BedTool(gtf_regions))

                print('*' * 72)
                line = '\t'.join(map(str, [species, release] + [result[o] for o in self.order]))
                print(line)
                print('*' * 72)
                self.current_results.append(line)

        # Confirm that reference file exists!
        self.assertTrue(isfile(self.ref_results_file))

        with open(self.ref_results_file, 'r') as ref_file:
            for reference_line in ref_file:
                self.reference_results.append(reference_line.strip())

        # Confirm that self.reference_results and self.current_results have
        # identical content:
        self.assertEqual(set(self.reference_results), set(self.current_results))

        self.success = True


if __name__ == '__main__':
    unittest.main()
