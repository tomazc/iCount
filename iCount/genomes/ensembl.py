"""
Ensembl API
-----------

Functions to query and download data from the `Ensembl`_ FTP site.

.. _Ensembl:
    http://www.ensembl.org/index.html

"""

import os
import re
import ftplib
import gzip
import shutil
import logging
import tempfile
import subprocess

import iCount

BASE_URL = 'ftp.ensembl.org'

MIN_RELEASE_SUPPORTED = 59
MAX_RELEASE_SUPPORTED = 84

LOGGER = logging.getLogger(__name__)


def _docstring_parameter(*sub):
    """
    To be used as decorator for passing arguments for docstring formatting.
    """
    def dec(obj):
        obj.__doc__ = obj.__doc__.format(*sub)
        return obj
    return dec


@_docstring_parameter(BASE_URL)
def get_ftp_instance():
    """
    Get ftplib.FTP object that is connected to {0}

    Returns
    -------
    ftplib.FTP
        FTP object connected to {0}
    """
    ftp = ftplib.FTP(BASE_URL)
    ftp.login()
    return ftp


@_docstring_parameter(MIN_RELEASE_SUPPORTED, MAX_RELEASE_SUPPORTED)
def releases():
    """
    Get list of available ENSEMBL releases.

    Only allows ENSEMBL releases from {0} - {1}.

    Returns
    -------
    list
        List of available releases

    """
    ftp = get_ftp_instance()
    # set current working directory
    ftp.cwd('pub')

    fetch = [int(item.strip('release-')) for item in ftp.nlst() if re.match(r'release-\d+', item)]
    ftp.quit()
    out = [i for i in fetch if i >= MIN_RELEASE_SUPPORTED and i <= MAX_RELEASE_SUPPORTED]
    out = sorted(out, reverse=True)

    LOGGER.info('There are %d releases available: %s', len(out), ','.join(map(str, (out))))
    return out


@_docstring_parameter(MIN_RELEASE_SUPPORTED, MAX_RELEASE_SUPPORTED)
def species(release):
    """
    Get list of species for given release.

    Parameters
    ----------
    release : str
        The release number (can be str or int). Only ENSEMBL releases
        from {0} - {1} are available.

    Returns
    -------
    list
        List of species.

    """
    iCount.log_inputs(LOGGER, level=logging.INFO)

    ftp = get_ftp_instance()
    ftp.cwd('pub/' + 'release-' + str(release) + '/fasta/')
    spec_list = sorted([item for item in ftp.nlst()])
    ftp.quit()
    LOGGER.info('There are %d species available: %s',
                len(spec_list), ','.join(map(str, spec_list)))
    return spec_list


@_docstring_parameter(MIN_RELEASE_SUPPORTED, MAX_RELEASE_SUPPORTED)
def annotation(release, species, target_dir=None, target_fname=None):
    """
    Download annotation in GTF file fromat.

    Parameters
    ----------
    release : int
        Release number. Only ENSEMBL releases from {0} - {1} are available.
    species : str
        Species latin name.
    target_dir : str
        Download to this directory (if not given, current working directory).
    target_fname : str
        Annotation filename (must have .gz file extension). If not given,
        species.release.gtf.gz is used.

    Returns
    -------
    str
        Downloaded annotation filename.

    TODO: target_dir & target_fname into one parameter? But the default name is
    a quite useful feature...
    """
    iCount.log_inputs(LOGGER, level=logging.INFO)

    if not target_dir:
        target_dir = os.getcwd()
    if not os.path.isdir(target_dir):
        raise ValueError('Directory "{}" does not exist'.format(target_dir))
    # Process filename
    if not target_fname:
        target_fname = '{}.{}.gtf.gz'.format(species, release)

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
            break

    if not annotation_file:
        LOGGER.info('No GTF file found for species %s, release %s', species, str(release))
        return None

    # Download to file on disk
    saveas_fname = os.path.join(target_dir, target_fname)
    LOGGER.info('Downloading annotation to: %s', saveas_fname)
    with open(saveas_fname, 'wb') as fhandle:
        ftp.retrbinary('RETR ' + annotation_file, fhandle.write)
    ftp.quit()

    LOGGER.info('Done.')
    return os.path.abspath(saveas_fname)


def chrom_length(fasta_in):
    """
    Compute chromosome lengths to a file by using samtools faidx.

    More about the .fai file format can be found here:
    http://www.htslib.org/doc/faidx.html

    Parameters
    ----------
    fasta_in : str
        Path to genome fasta file (can be a .gz file)

    Returns
    -------
    str
        Absoulute path to output file

    """
    iCount.log_inputs(LOGGER, level=logging.INFO)

    temp = iCount.files.decompress_to_tempfile(fasta_in)
    command = ['samtools', 'faidx', temp]
    subprocess.check_call(command)

    # This command makes fai file. Move & rename this file to fasta_in.fai:
    fai_file_temp = temp + '.fai'
    fai_file = fasta_in + '.fai'
    subprocess.check_call(['mv', fai_file_temp, fai_file])
    LOGGER.info('Fai file just saved to : %s', os.path.abspath(fai_file))
    return os.path.abspath(fai_file)


@_docstring_parameter(MIN_RELEASE_SUPPORTED, MAX_RELEASE_SUPPORTED)
def sequence(release, species, target_dir=None, target_fname=None, tempdir=None, chromosomes=None):
    """
    Downloads genome file in FASTA fromat.

    Several steps are performed:

        * querry for list off all FASTA files for given release and species
        * filter this list to get only whole chromosome files
        * if chromosomes paramter is given, take only specified chromosomes
        * sort list of these files to have the correct order
        * download each file and write it to target_fname

    Parameters
    ----------
    release : int
        Release number. Only ENSEMBL releases from {0} - {1} are available.
    species : str
        Species latin name.
    target_dir : str
        Download to this directory (if not given, current working directory).
    target_fname : str
        Annotation filename (must have .gz file extension). If not given,
        species.release.gtf.gz is used.
    tempdir : str
        Temporary directory with intermediate results.
    chromosomes : list_str
        If given, don't download the whole genome, but juts the given
        cromosomes. Chromosomes can be given as strings aor integers

    Returns
    -------
    str
        Downloaded genome/sequnce filename.

    TODO: target_dir & target_fname into one parameter? But the default name is
    a quite useful feature...
    """
    iCount.log_inputs(LOGGER, level=logging.INFO)

    if not target_dir:
        target_dir = os.getcwd()
    if not os.path.isdir(target_dir):
        raise ValueError('Directory "{}" does not exist'.format(target_dir))

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
        if not target_fname:
            target_fname = '{}.{}.chr{}.fa.gz'.format(
                species, release, '_'.join(map(str, chromosomes)))

    if not target_fname:
        target_fname = '{}.{}.fa.gz'.format(species, release)
    target_path = os.path.abspath(os.path.join(target_dir, target_fname))

    LOGGER.info('Downloading FASTA file into: %s', target_path)
    tempdir = tempfile.mkdtemp(dir=tempdir)
    for fname in filtered_files:
        LOGGER.debug('Downloading file: %s', fname)
        # Download all files to tempdir:
        temp_file = os.path.join(tempdir, fname)
        with open(temp_file, 'wb') as fhandle:
            ftp.retrbinary('RETR ' + fname, fhandle.write)

        # Write the content into final file (target_fname)
        with gzip.open(temp_file, 'rb') as f_in, gzip.open(target_path, 'ab') as f_out:
            shutil.copyfileobj(f_in, f_out)

    # Clean up: delete temporary files:
    ftp.quit()
    shutil.rmtree(tempdir)

    # Compute chromosome lengths:
    chrom_length(target_path)

    LOGGER.info('Done.')
    return target_path
