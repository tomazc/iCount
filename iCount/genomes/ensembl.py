"""
Functions to querry and download data from ENSEMBL ftp site
"""

import os
import re
import uuid
import ftplib
import gzip
import shutil


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
    # TODO Jure: list should be sorted by decreasing version (latest first)
    # TODO Jure: limit releases from ensembl 59 (including) onwards
    # TODO Jure: limit upwards as well, use MIN_RELEASE_SUPPORTED and MAX_RELEASE_SUPPORTED
    ftp = get_ftp_instance()
    # set current working directory
    ftp.cwd('pub')
    # ftp.nlst() returns a list of all files/folders in cwd
    return [item.strip('release-') for item in ftp.nlst()
            if re.match(r'release-\d+', item)]


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
    if species_list:
        return species_list
    else:
        raise ValueError('No specias found for release {}'.format(release))


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

    # TODO Jure: return None if file not found on FTP

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

    # Process filename
    if not target_fname:
        target_fname = '{}.{}.gtf.gz'.format(species, release)

    saveas_fname = os.path.join(target_dir, target_fname)

    # Download to file on disk
    with open(saveas_fname, 'wb') as fhandle:
        ftp.retrbinary('RETR ' + annotation_file, fhandle.write)
    return saveas_fname


def download_sequence(release, species, target_dir=None, target_fname=None,
                      tempdir=None):
    """
    Downloads whole genome file for given release and species

    Several steps are performed:

        * querry for list off all FASTA files for given release and specias
        * filter this list to get only whole chromosome files
        * sort list of these files to have the correct order
        * download each file and write it to target_fname

    :param str/int release: the release number
    :param str species: the species latin name
    :param str/path target_dir: download location
    :param str target_fname: desired name of downloaded file (must have .gz file extension)
    :param str tempdir: location of temp. folder and files

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
    target_path = os.path.join(target_dir, target_fname)

    ftp = get_ftp_instance()
    server_dir = '/pub/release-{}/fasta/{}/dna'.format(release, species)
    ftp.cwd(server_dir)
    all_fasta_files = ftp.nlst()

    filterd_files = []

    # Recognize all "whole-chromosme" files
    regex = r'{}\.[\w\d]+\.(\d+\.)*dna\.chromosome\.' \
            r'[\dXYMT]+\.fa\.gz'.format(species.capitalize())
    for fasta_file in all_fasta_files:
        if re.match(regex, fasta_file):
            filterd_files.append(fasta_file)

    # Make sure there is temporary directory to store temporary files
    # TODO Jure: there could be other stuff in the specified tempdir, always create a temporary subfolder in the tempdir
    if not tempdir:
        # gnerate ranodm folder name to avoid clashes with existing folders:
        tempdir = os.path.join(os.getcwd(), str(uuid.uuid4()))
        os.mkdir(tempdir)

    def sorting_func(name):
        """Helper function for sorting files named by chromosome"""
        key = name.split('.')[-3]
        if key.isdigit():
            return int(key)
        elif key == 'MT':
            return 'ZZ'  # this makes mitohondrial "chromosome" the last one
        return key
    # Sort the filtered_files in correct order of chromosomes:
    filterd_files.sort(key=sorting_func)

    for fname in filterd_files:
        # Download all files to tempdir:
        tempfile = os.path.join(tempdir, fname)
        with open(tempfile, 'wb') as fhandle:
            ftp.retrbinary('RETR ' + fname, fhandle.write)

        # Write the content into final file (target_fname)
        with gzip.open(tempfile, 'rb') as f_in, \
             gzip.open(target_path, 'ab') as f_out:
            shutil.copyfileobj(f_in, f_out)

    # Clean up: delete temporary files:
    shutil.rmtree(tempdir)


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
