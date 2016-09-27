import os
import unittest
import ftplib
import tempfile

from iCount.genomes import ensembl


class TestEnsembl(unittest.TestCase):

    def test_get_ftp_instance(self):
        ftp_instance = ensembl.get_ftp_instance()

        self.assertIsInstance(ftp_instance, ftplib.FTP)
        self.assertIsInstance(ftp_instance.pwd(), str)

        ftp_instance.quit()

    def test_get_release_list(self):
        releases = ensembl.get_release_list()

        self.assertIsInstance(releases, list)
        self.assertTrue(len(releases) > 0)
        self.assertIsInstance(releases[0], int)
        self.assertEqual(min(releases), 59)
        self.assertTrue(max(releases) > 83)

    def test_get_species_list(self):
        species = ensembl.get_species_list(84)

        self.assertIsInstance(species, list)
        self.assertTrue(len(species) > 0)
        self.assertIsInstance(species[0], str)
        self.assertTrue('homo_sapiens' in species)


class TestEnsemblDownload(unittest.TestCase):

    def setUp(self):
        self.tempdir = tempfile.mkdtemp()
        self.gtf = 'test_file.gtf'
        self.fasta = 'test_file.fa.gz'

    def test_download_annotation(self):
        ann_file = ensembl.download_annotation(
            '84', 'homo_sapiens', target_dir=self.tempdir, target_fname=self.gtf)

        self.assertTrue(os.path.isfile(os.path.join(self.tempdir, self.gtf)))

    def test_download_sequence(self):
        fasta_file = ensembl.download_sequence(
            '84', 'homo_sapiens', target_dir=self.tempdir, chromosomes=['MT'])

        self.assertTrue(os.path.isfile(
            os.path.join(self.tempdir, 'homo_sapiens.84.chrMT.fa.gz')))
        # Confirm that chrom_length file was created!
        self.assertTrue(os.path.isfile(
            os.path.join(self.tempdir, 'homo_sapiens.84.chrMT.fa.gz.fai')))

    def tearDown(self):
        files = os.listdir(self.tempdir)
        for file_ in files:
            os.remove(os.path.join(self.tempdir, file_))

        os.rmdir(self.tempdir)

if __name__ == '__main__':
    unittest.main()
