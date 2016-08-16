"""
This script tests if results of iCount.genomes.segments.get_genes function
are consistent with the reference results. If reference results are not
present, they are created.

If gtf file for given release and species is not found in `self.gtf_files_dir`
it is downloaded from ftp.ensembl.org by using function
iCount.genomes.ensembl.download_annotation. This may significantly extend
the execution time of this script.

By modifying the variables `self.releases_list` and `self.species_list` user can
determine which gtf files to put under test.
"""

import os
import unittest

from uuid import uuid4

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

        self.current_results = []
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
                gtf_in = os.path.join(self.gtf_files_dir, '{}.{}.gtf.gz'.format(species, release))
                # If gtf_in file is not preswent, download it:
                if not os.path.isfile(gtf_in):
                    _ = ensembl.download_annotation(
                        release, species, target_dir=self.gtf_files_dir)

                gtf_out = os.path.join(
                    self.gtf_files_dir, '{}.{}.genes.gtf.gz'.format(species, str(release)))

                gene_segments = segment.get_genes(gtf_in, gtf_out, attribute='gene_id')
                line = 'species: {}, release: {}, segments: {}'.format(
                    species, release, len(gene_segments))
                print(line)
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
