# pylint: disable=missing-docstring, protected-access

import ftplib
import os
import unittest
from unittest import mock
import warnings

from iCount import genomes
from iCount.genomes import ensembl
from iCount.tests.utils import get_temp_file_name, get_temp_dir


class TestEnsemblUtils(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", ResourceWarning)

    def test_get_ftp_instance(self):
        ftp_instance = genomes.get_ftp_instance(ensembl.BASE_URL)

        self.assertIsInstance(ftp_instance, ftplib.FTP)
        self.assertIsInstance(ftp_instance.pwd(), str)

        ftp_instance.quit()

    @mock.patch('iCount.genomes.ftplib')
    def test_no_connection(self, ftplib_mock):
        ftplib_mock.FTP = mock.MagicMock(side_effect=Exception())
        message = "Problems connecting to ENSEMBL FTP server."
        with self.assertRaisesRegex(Exception, message):
            genomes.get_ftp_instance(ensembl.BASE_URL)


class TestSpecies(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", ResourceWarning)

    def test_bad_source(self):
        message = r'Source .* is not supported.'
        with self.assertRaisesRegex(ValueError, message):
            genomes.species(source='invalid_source')

    def test_ensembl_1(self):
        species = genomes.species(source='ensembl')

        self.assertIsInstance(species, list)
        self.assertTrue(len(species) > 0)
        self.assertTrue('homo_sapiens' in species)

    def test_ensembl_2(self):
        species = genomes.species(source='ensembl', release=84)

        self.assertIsInstance(species, list)
        self.assertTrue(len(species) > 0)
        self.assertTrue('homo_sapiens' in species)

    def test_gencode(self):
        species = genomes.species()

        self.assertIsInstance(species, list)
        self.assertTrue(len(species) > 0)
        self.assertIsInstance(species[0], str)
        self.assertTrue('human' in species)
        self.assertTrue('mouse' in species)


class TestReleases(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", ResourceWarning)

    def test_bad_source(self):
        message = r'Source .* is not supported.'
        with self.assertRaisesRegex(ValueError, message):
            genomes.releases(source='invalid_source')

    def test_ensembl(self):
        releases = genomes.releases(source='ensembl')

        self.assertIsInstance(releases, list)
        self.assertTrue(len(releases) > 0)
        self.assertTrue(84 in releases)

    def test_gencode(self):
        releases = genomes.releases(source='gencode', species='human')

        self.assertIsInstance(releases, list)
        self.assertTrue(len(releases) > 0)
        self.assertIsInstance(releases[0], str)
        self.assertTrue('27' in releases)


class TestEnsemblAnnotation(unittest.TestCase):

    def setUp(self):
        self.tempdir = get_temp_dir()
        self.tmpfile = get_temp_file_name(extension='.gtf.gz')
        warnings.simplefilter("ignore", ResourceWarning)

    def test_annotation(self):
        genomes.annotation('homo_sapiens', release=84, out_dir=self.tempdir,
                           annotation=self.tmpfile, source='ensembl')
        self.assertTrue(os.path.isfile(os.path.join(self.tempdir, self.tmpfile)))

    def test_annotation_invalid_species(self):
        message = r'Invalid species name.'
        with self.assertRaisesRegex(ValueError, message):
            genomes.annotation('invalid_species_name', 84, source='ensembl')

    def test_annotation_invalid_release(self):
        message = r"Release should be a number between \d+ and \d+"
        with self.assertRaisesRegex(ValueError, message):
            genomes.annotation('homo_sapiens', 1, source='ensembl')

    def test_annotation_invalid_out_dir(self):
        message = r"Directory .* does not exist."
        with self.assertRaisesRegex(ValueError, message):
            genomes.annotation('homo_sapiens', release=84, out_dir='/out/dir/that/does/not/exist',
                               source='ensembl')

    def tearDown(self):
        files = os.listdir(self.tempdir)
        for file_ in files:
            os.remove(os.path.join(self.tempdir, file_))

        os.rmdir(self.tempdir)


class TestGencodeAnnotation(unittest.TestCase):

    def setUp(self):
        self.tempdir = get_temp_dir()
        self.tmpfile = get_temp_file_name(extension='gtf.gz')
        warnings.simplefilter("ignore", ResourceWarning)

    def test_bad_source(self):
        message = r'Source .* is not supported.'
        with self.assertRaisesRegex(ValueError, message):
            genomes.annotation('human', '27', source='invalid_source')

    def test_annotation(self):
        genomes.annotation('mouse', release='M15', out_dir=self.tempdir, annotation=self.tmpfile,
                           source='gencode')
        self.assertTrue(os.path.isfile(os.path.join(self.tempdir, self.tmpfile)))

    def test_annotation_invalid_species(self):
        message = r'Invalid species name.'
        with self.assertRaisesRegex(ValueError, message):
            genomes.annotation(species='invalid_species_name', release='26', source='gencode')

    def test_annotation_invalid_release(self):
        message = r"Invalid release number."
        with self.assertRaisesRegex(ValueError, message):
            genomes.annotation('human', release='42000', source='gencode')

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
        genomes.genome(
            'homo_sapiens', release=84, out_dir=self.tempdir, chromosomes=['MT'], source='ensembl')

        self.assertTrue(os.path.isfile(
            os.path.join(self.tempdir, 'homo_sapiens.84.chrMT.fa.gz')))
        # Confirm that chrom_length file was created!
        self.assertTrue(os.path.isfile(
            os.path.join(self.tempdir, 'homo_sapiens.84.chrMT.fa.gz.fai')))

    def test_genome_invalid_species(self):
        message = r'Invalid species name.'
        with self.assertRaisesRegex(ValueError, message):
            genomes.genome('invalid_species_name', 84, source='ensembl')

    def test_genome_invalid_release(self):
        message = r"Release should be a number between \d+ and \d+"
        with self.assertRaisesRegex(ValueError, message):
            genomes.genome('homo_sapiens', release=1, source='ensembl')

    def test_invalid_chromosome(self):
        message = r"Could not find chromosome .*"
        with self.assertRaisesRegex(ValueError, message):
            genomes.genome('homo_sapiens', 84, chromosomes=['ZZZ'], source='ensembl')

    def test_invalid_out_dir(self):
        message = r'Directory ".*" does not exist.'
        with self.assertRaisesRegex(ValueError, message):
            genomes.genome('homo_sapiens', 84, genome='genome.fa.gz', out_dir='/does/not/exist',
                           source='ensembl')

    def test_absolute_path(self):
        message = r'Directory ".*" does not exist.'
        with self.assertRaisesRegex(ValueError, message):
            genomes.genome('homo_sapiens', 84, genome='/invalid/path/genome.fa.gz',
                           out_dir='/invalid/path/will/raise/error', source='ensembl')

    def tearDown(self):
        files = os.listdir(self.tempdir)
        for file_ in files:
            os.remove(os.path.join(self.tempdir, file_))

        os.rmdir(self.tempdir)


class TestGencodeGenome(unittest.TestCase):

    def setUp(self):
        self.tempdir = get_temp_dir()
        warnings.simplefilter("ignore", ResourceWarning)

    def test_bad_source(self):
        message = r'Source .* is not supported.'
        with self.assertRaisesRegex(ValueError, message):
            genomes.genome('human', '27', source='invalid_source')

    @unittest.skip("This file is too large to download it in unit test.")
    def test_genome(self):
        genomes.genome('human', release='26', out_dir=self.tempdir, source='gencode')

        self.assertTrue(os.path.isfile(os.path.join(self.tempdir, 'human.26.fa.gz')))
        # Confirm that chrom_length file was created!
        self.assertTrue(os.path.isfile(os.path.join(self.tempdir, 'human.26.fa.gz.fai')))

    def test_genome_invalid_species(self):
        message = r'Invalid species name.'
        with self.assertRaisesRegex(ValueError, message):
            genomes.genome(species='invalid_species_name', release='26', source='gencode')

    def test_genome_invalid_release(self):
        message = r'Invalid release number.'
        with self.assertRaisesRegex(ValueError, message):
            genomes.genome('human', release='1000', source='gencode')

    def tearDown(self):
        files = os.listdir(self.tempdir)
        for file_ in files:
            os.remove(os.path.join(self.tempdir, file_))

        os.rmdir(self.tempdir)


if __name__ == '__main__':
    unittest.main()
