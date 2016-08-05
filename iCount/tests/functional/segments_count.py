"""
This script tests if results of iCount.genomes.segments.get_regions function
are consistent with the reference results. If reference results are not
present, they are created.

If gtf/fasta/genome/get_regions.gtf files for given release and species are
not found in `self.gtf_files_dir` they are downloaded/calculated. This may
significantly extend the execution time of this script.

By modifying the variables `self.releases_list` and `self.species_list` user can
determine which gtf files to put under test.
"""

import os
from uuid import uuid4
import unittest

import pybedtools

from collections import Counter

from iCount.genomes import ensembl, segment
from __init__ import functional_test_folder


class TestSegmentation(unittest.TestCase):

    def setUp(self):
        self.ref_results_file = __file__[:-3] + '.out.txt'

        # Folder where all gtf files are stored:
        self.gtf_files_dir = os.path.join(functional_test_folder, 'gtf_files/')
        if not os.path.exists(self.gtf_files_dir):
            os.makedirs(self.gtf_files_dir)

        self.releases_list = ensembl.get_release_list()
        self.species_list = ['homo_sapiens', 'mus_musculus']

        self.order = ['gene', 'transcript', 'CDS', 'UTR3', 'UTR5', 'intron',
                      'stop_codon', 'ncRNA', 'intergenic']
        header = '\t'.join(map(str, ['species', 'release'] + [o for o in self.order]))

        self.current_results = [header]
        self.reference_results = []
        self.success = False

    def tearDown(self):
        if not self.success and self.current_results:
            if os.path.isfile(self.ref_results_file):
                self.ref_results_file = self.ref_results_file + str(uuid4())
            with open(self.ref_results_file, 'w') as ref_file:
                for result in self.current_results:
                    ref_file.write(result + '\n')

    def test_reference_stats(self):
        for species in self.species_list:
            for release in self.releases_list:
                gtf_regions = os.path.join(
                    self.gtf_files_dir,
                    '{}.{}.get_regions.gtf.gz'.format(species, release))

                if not os.path.isfile(gtf_regions):
                    # If gtf_regions file is not preswent, make it from
                    # gtf_base and chrom_length file
                    gtf_base = os.path.join(self.gtf_files_dir, '{}.{}.gtf.gz'.format(species, release))
                    if not os.path.isfile(gtf_base):
                        _ = ensembl.download_annotation(release, species, target_dir=self.gtf_files_dir)

                    chrom_length = os.path.join(
                        self.gtf_files_dir, '{}.{}.fa.gz.chrom_length.txt'.format(species, release))
                    if not os.path.isfile(chrom_length):
                        fasta = os.path.join(self.gtf_files_dir, '{}.{}.fa.gz'.format(species, release))
                        if not os.path.isfile(fasta):
                            fasta = ensembl.download_sequence(release, species, target_dir=self.gtf_files_dir)
                        # Now we have all the fasta file to produce chrom_length
                        chrom_length = ensembl.chrom_length(fasta)

                    # Now we surely have gtf_base and chrom_length to produce gtf_regions:
                    gtf_regions = segment.get_regions(gtf_base, gtf_regions, chrom_length)
                    # When multiple-cpu usage is supported in get_regions, one can
                    # speed up the test by using multiple cpus like this:
                    # gtf_regions = segment.get_regions(gtf_base, gtf_regions, chrom_length, cores=6)

                result = Counter(i[2] for i in pybedtools.BedTool(gtf_regions))

                print('*' * 72)
                line = '\t'.join(map(str, [species, release] + [result[o] for o in self.order]))
                print(line)
                print('*' * 72)
                self.current_results.append(line)

        # Confirm that reference file exists!
        self.assertTrue(os.path.isfile(self.ref_results_file))

        with open(self.ref_results_file, 'r') as ref_file:
            for reference_line in ref_file:
                self.reference_results.append(reference_line.strip())

        # Confirm that self.reference_results and self.current_results have
        # identical content:
        self.assertEqual(set(self.reference_results), set(self.current_results))

        self.success = True


if __name__ == '__main__':
    unittest.main()
