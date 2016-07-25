"""
Functions to querry and download data from ENSEMBL ftp site
"""

import os
import re
import ftplib
import gzip
import shutil
import tempfile


BASE_URL = 'ftp.ensembl.org'

MIN_RELEASE_SUPPORTED = 59
MAX_RELEASE_SUPPORTED = 84


def get_ftp_instance():
    """
    Get ftplib.FTP object that is connected to BASE_URL

    :return: FTP object connected to BASE_URL
    :rtype: ftplib.FTP
    """
    ftp = ftplib.FTP(BASE_URL)
    ftp.login()
    return ftp


def get_release_list():
    """
    Get list of available releases

    :return: list of available releases - elements are numbers of type str
    :rtype: list
    """
    ftp = get_ftp_instance()
    # set current working directory
    ftp.cwd('pub')

    out = [item.strip('release-') for item in ftp.nlst() if
           re.match(r'release-\d+', item) and
           int(item.strip('release-')) >= MIN_RELEASE_SUPPORTED and
           int(item.strip('release-')) <= MAX_RELEASE_SUPPORTED]

    ftp.quit()
    return sorted(out, reverse=True)


def get_species_list(release):
    """
    Get list of species for given release

    :param str release: the release number of type str or int
    :return: list of species - elements are latin_names
    :rtype: list
    :raises ValueError: if no species are found for given release
    """
    ftp = get_ftp_instance()
    ftp.cwd('pub/' + 'release-' + str(release) + '/fasta/')
    species_list = [item for item in ftp.nlst()]

    ftp.quit()
    if species_list:
        return sorted(species_list)


def download_annotation(release, species, target_dir=None, target_fname=None):
    """
    Downloads annotation file for given release and species

    :param str/int release: the release number
    :param str species: the species latin name
    :param str/path target_dir: download location
    :param str target_fname: desired name of downloaded file (must have .gz file extension)

    :return: filename of downloaded file
    :rtype: str
    :raises ValueError: if provided target_dir does not exist
    """

    if not target_dir:
        target_dir = os.getcwd()
    if not os.path.isdir(target_dir):
        raise ValueError('Directory "{}" does not exist'.format(target_dir))

    ftp = get_ftp_instance()
    server_dir = '/pub/release-{}/gtf/{}/'.format(release, species)
    ftp.cwd(server_dir)
    server_files = ftp.nlst()

    # Get the desired annotation file form list of all server files:
    annotation_file = ''
    regex = r'{}\.[\w\d]+\.\d+\.gtf\.gz'.format(species.capitalize())
    for file_ in server_files:
        if re.match(regex, file_):
            annotation_file = file_

    if not annotation_file:
        return None

    # Process filename
    if not target_fname:
        target_fname = '{}.{}.gtf.gz'.format(species, release)

    saveas_fname = os.path.join(target_dir, target_fname)

    # Download to file on disk
    with open(saveas_fname, 'wb') as fhandle:
        ftp.retrbinary('RETR ' + annotation_file, fhandle.write)

    ftp.quit()
    return saveas_fname


def download_sequence(release, species, target_dir=None, target_fname=None,
                      tempdir=None, chromosomes=[]):
    """
    Downloads whole genome file for given release and species

    Several steps are performed:

        * querry for list off all FASTA files for given release and specias
        * filter this list to get only whole chromosome files
        * if chromosomes paramter is given, take only specified chromosomes
        * sort list of these files to have the correct order
        * download each file and write it to target_fname

    :param str/int release: the release number
    :param str species: the species latin name
    :param str/path target_dir: download location
    :param str target_fname: desired name of downloaded file (must have .gz file extension)
    :param str tempdir: location of temp. folder and files
    :param list chromosomes: a subset of chromosomes to download

    :return: filename of downloaded file
    :rtype: str
    :raises ValueError: if provided target_dir does not exist
    """

    if not target_dir:
        target_dir = os.getcwd()
    if not os.path.isdir(target_dir):
        raise ValueError('Directory "{}" does not exist'.format(target_dir))
    if not target_fname:
        target_fname = '{}.{}.fa.gz'.format(species, release)

    ftp = get_ftp_instance()
    server_dir = '/pub/release-{}/fasta/{}/dna'.format(release, species)
    ftp.cwd(server_dir)
    all_fasta_files = ftp.nlst()

    filtered_files = []

    # Recognize all "whole-chromosme" files
    regex = r'{}\.[\w\d]+\.(\d+\.)*dna\.chromosome\.' \
            r'[\dXYMT]+\.fa\.gz'.format(species.capitalize())
    for fasta_file in all_fasta_files:
        if re.match(regex, fasta_file):
            filtered_files.append(fasta_file)

    def sorting_func(name):
        """Helper function for sorting files named by chromosome"""
        key = name.split('.')[-3]
        if key == 'MT':
            return 'ZZ'  # this makes mitohondrial "chromosome" the last one
        return key
    # Sort the filtered_files in correct order of chromosomes:
    filtered_files.sort(key=sorting_func)

    # If parameter chromosomes is given, take only a subset of all chromosomes
    if chromosomes:
        subset_list = []
        for file_ in filtered_files:
            for chromosome in chromosomes:
                regex = r'.*chromosome\.{}\.fa\.gz'.format(str(chromosome))
                if re.match(regex, file_):
                    subset_list.append(file_)
                else:
                    pass
        filtered_files = subset_list
        target_fname = target_fname.rstrip('.fa.gz') + '.chr' + \
            '_'.join(map(str, chromosomes)) + '.fa.gz'

    tempdir = tempfile.mkdtemp(dir=tempdir)
    target_path = os.path.join(target_dir, target_fname)

    for fname in filtered_files:
        # Download all files to tempdir:
        temp_file = os.path.join(tempdir, fname)
        with open(temp_file, 'wb') as fhandle:
            ftp.retrbinary('RETR ' + fname, fhandle.write)

        # Write the content into final file (target_fname)
        with gzip.open(temp_file, 'rb') as f_in, \
             gzip.open(target_path, 'ab') as f_out:
            shutil.copyfileobj(f_in, f_out)

    # Clean up: delete temporary files:
    ftp.quit()
    shutil.rmtree(tempdir)

    return target_path


# #############################################################

if __name__ == '__main__':

    releases = get_release_list()
    release84 = releases[-2]
    species = get_species_list(releases[-2])

    # Choose random species
    species = species[42]

    # download_annotation(release84, species, '/home/jure/Desktop')

    # Executing this command might take a looong time!
    # download_genome(release84, species, '/home/jure/Desktop')
