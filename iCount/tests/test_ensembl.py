# pylint: disable=missing-docstring, protected-access

import os
import ftplib
import warnings
import unittest
from unittest import mock

from iCount.genomes import ensembl
from iCount.tests.utils import get_temp_file_name, get_temp_dir


class TestEnsembl(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", ResourceWarning)

    def test_get_ftp_instance(self):
        ftp_instance = ensembl.get_ftp_instance()

        self.assertIsInstance(ftp_instance, ftplib.FTP)
        self.assertIsInstance(ftp_instance.pwd(), str)

        ftp_instance.quit()

    @mock.patch('iCount.genomes.ensembl.ftplib')
    def test_no_connection(self, ftplib_mock):
        ftplib_mock.FTP = mock.MagicMock(side_effect=Exception())
        message = "Problems connecting to ENSEMBL FTP server."
        with self.assertRaisesRegex(Exception, message):
            ensembl.get_ftp_instance()

    def test_releases(self):
        releases = ensembl.releases()

        self.assertIsInstance(releases, list)
        self.assertTrue(len(releases) > 0)
        self.assertIsInstance(releases[0], int)
        self.assertEqual(min(releases), 59)
        self.assertTrue(max(releases) > 83)

    def test_species(self):
        species = ensembl.species()

        self.assertIsInstance(species, list)
        self.assertTrue(len(species) > 0)
        self.assertIsInstance(species[0], str)
        self.assertTrue('homo_sapiens' in species)

        species2 = ensembl.species(84)
        self.assertEqual(species2, species)

    def test_species_bad_release_number(self):
        message = r"Release should be a number between \d+ and \d+"
        with self.assertRaisesRegex(ValueError, message):
            ensembl.species(1)


class TestEnsemblAnnotation(unittest.TestCase):

    def setUp(self):
        self.tempdir = get_temp_dir()
        self.tmpfile = get_temp_file_name(extension='.gtf.gz')
        warnings.simplefilter("ignore", ResourceWarning)

    def test_annotation(self):
        ensembl.annotation(
            'homo_sapiens', release=84, out_dir=self.tempdir, annotation=self.tmpfile)

        self.assertTrue(os.path.isfile(os.path.join(self.tempdir, self.tmpfile)))

    def test_annotation_invalid_species(self):
        message = r'Invalid species name.'
        with self.assertRaisesRegex(ValueError, message):
            ensembl.annotation('invalid_species_name')

    def test_annotation_invalid_release(self):
        message = r"Release should be a number between \d+ and \d+"
        with self.assertRaisesRegex(ValueError, message):
            ensembl.annotation('homo_sapiens', release=1)

    def tearDown(self):
        files = os.listdir(self.tempdir)
        for file_ in files:
            os.remove(os.path.join(self.tempdir, file_))

        os.rmdir(self.tempdir)


class TestEnsemblGenome(unittest.TestCase):

    def setUp(self):
        self.tempdir = get_temp_dir()
        warnings.simplefilter("ignore", ResourceWarning)

    def test_genome(self):
        ensembl.genome(
            'homo_sapiens', release=84, out_dir=self.tempdir, chromosomes=['MT'])

        self.assertTrue(os.path.isfile(
            os.path.join(self.tempdir, 'homo_sapiens.84.chrMT.fa.gz')))
        # Confirm that chrom_length file was created!
        self.assertTrue(os.path.isfile(
            os.path.join(self.tempdir, 'homo_sapiens.84.chrMT.fa.gz.fai')))

    def test_genome_invalid_species(self):
        message = r'Invalid species name.'
        with self.assertRaisesRegex(ValueError, message):
            ensembl.genome('invalid_species_name')

    def test_genome_invalid_release(self):
        message = r"Release should be a number between \d+ and \d+"
        with self.assertRaisesRegex(ValueError, message):
            ensembl.genome('homo_sapiens', release=1)

    def tearDown(self):
        files = os.listdir(self.tempdir)
        for file_ in files:
            os.remove(os.path.join(self.tempdir, file_))

        os.rmdir(self.tempdir)


if __name__ == '__main__':
    unittest.main()
