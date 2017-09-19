""".. Line to protect from pydocstyle D205, D400.

Ensembl API
-----------

Functions to query and download data from the `Ensembl`_ FTP site.

.. _Ensembl:
    http://www.ensembl.org/index.html

"""

import gzip
import logging
import os
import re
import shutil

import pysam

import iCount

BASE_URL = 'ftp.ensembl.org'

MIN_RELEASE_SUPPORTED = 59
MAX_RELEASE_SUPPORTED = 88

LOGGER = logging.getLogger(__name__)


def _docstring_parameter(*sub):
    """Pass arguments to function docstrings (to be used as decorator)."""
    def dec(obj):  # pylint: disable=missing-docstring
        obj.__doc__ = obj.__doc__.format(*sub)
        return obj
    return dec


@_docstring_parameter(MIN_RELEASE_SUPPORTED, MAX_RELEASE_SUPPORTED)
def releases():
    """
    Get list of available ENSEMBL releases.

    Only allows ENSEMBL releases {0}-{1}.

    Returns
    -------
    list
        List of available releases

    """
    ftp = iCount.genomes.get_ftp_instance(BASE_URL)
    # set current working directory
    ftp.cwd('pub')

    fetch = [int(item.strip('release-')) for item in ftp.nlst() if re.match(r'release-\d+', item)]
    ftp.quit()
    out = [i for i in fetch if i >= MIN_RELEASE_SUPPORTED and i <= MAX_RELEASE_SUPPORTED]
    out = sorted(out, reverse=True)

    LOGGER.info('There are %d ENSEMBL releases available: %s', len(out), ','.join(map(str, (out))))
    return out


@_docstring_parameter(MIN_RELEASE_SUPPORTED, MAX_RELEASE_SUPPORTED)
def species(release=MAX_RELEASE_SUPPORTED):
    """
    Get list of available species for given ENSEMBL release.

    Parameters
    ----------
    release : int
        Release number. Only ENSEMBL releases {0}-{1} are available.

    Returns
    -------
    list
        List of species.

    """
    if release is None:
        release = MAX_RELEASE_SUPPORTED
    iCount.log_inputs(LOGGER, level=logging.INFO)

    if not MIN_RELEASE_SUPPORTED <= release <= MAX_RELEASE_SUPPORTED:
        raise ValueError('Release should be a number between {} and {}.'.format(
            MIN_RELEASE_SUPPORTED, MAX_RELEASE_SUPPORTED))

    ftp = iCount.genomes.get_ftp_instance(BASE_URL)
    ftp.cwd('pub/' + 'release-' + str(release) + '/fasta/')
    spec_list = sorted([item for item in ftp.nlst()])
    ftp.quit()
    LOGGER.info('There are %d species available: %s',
                len(spec_list), ','.join(map(str, spec_list)))
    return spec_list


@_docstring_parameter(MIN_RELEASE_SUPPORTED, MAX_RELEASE_SUPPORTED)
# pylint: disable=redefined-outer-name
def annotation(species, release=MAX_RELEASE_SUPPORTED, out_dir=None, annotation=None):
    """
    Download ENSEMBL annotation for given release/species.

    Parameters
    ----------
    species : str
        Species latin name.
    release : int
        Release number. Only ENSEMBL releases {0}-{1} are available.
    out_dir : str
        Download to this directory (if not given, current working directory).
    annotation : str
        Annotation filename (must have .gz file extension). If not given,
        species.release.gtf.gz is used. If annotation is provided as absolute
        path, value of out_dir parameter is ignored and file is saved to given
        absolute path.

    Returns
    -------
    str
        Downloaded annotation filename.

    """
    iCount.log_inputs(LOGGER, level=logging.INFO)

    if not MIN_RELEASE_SUPPORTED <= release <= MAX_RELEASE_SUPPORTED:
        raise ValueError('Release should be a number between {} and {}.'.format(
            MIN_RELEASE_SUPPORTED, MAX_RELEASE_SUPPORTED))

    if species not in iCount.genomes.ensembl.species(release):
        raise ValueError('Invalid species name.')

    if annotation:
        assert annotation.endswith(('.gtf', '.gtf.gz'))
    else:
        annotation = '{}.{}.gtf.gz'.format(species, release)

    # If absolute path is given, ignore out_dir:
    if os.path.isabs(annotation):
        if out_dir:
            LOGGER.info('out_dir parameter has been changed, since absolute '
                        'path is provided by annotation parameter.')
        out_dir = os.path.dirname(annotation)
        annotation = os.path.basename(annotation)
    if not out_dir:
        out_dir = os.getcwd()
    if not os.path.isdir(out_dir):
        raise ValueError('Directory "{}" does not exist.'.format(out_dir))

    ftp = iCount.genomes.get_ftp_instance(BASE_URL)
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
    saveas_fname = os.path.join(out_dir, annotation)
    LOGGER.info('Downloading GTF to: %s', saveas_fname)
    with open(saveas_fname, 'wb') as fhandle:
        ftp.retrbinary('RETR ' + annotation_file, fhandle.write)
    ftp.quit()

    LOGGER.info('Done.')
    return os.path.abspath(saveas_fname)


def chrom_length(fasta_in):
    """
    Compute chromosome lengths of fasta file and store them into a file.

    More about the .fai file format can be found here:
    http://www.htslib.org/doc/faidx.html

    Parameters
    ----------
    fasta_in : str
        Path to genome FASTA file (can be .gz).

    Returns
    -------
    str
        Absolute path to output file.

    """
    iCount.log_inputs(LOGGER, level=logging.INFO)

    temp = iCount.files.decompress_to_tempfile(fasta_in)
    pysam.faidx(temp)  # pylint: disable=no-member

    fai_file = os.path.abspath(fasta_in + '.fai')
    shutil.move(temp + '.fai', fai_file)
    LOGGER.info('Fai file saved to : %s', fai_file)
    return fai_file


@_docstring_parameter(MIN_RELEASE_SUPPORTED, MAX_RELEASE_SUPPORTED)
# pylint: disable=redefined-outer-name
def genome(species, release=MAX_RELEASE_SUPPORTED, out_dir=None, genome=None,
           chromosomes=None):
    """
    Download ENSEMBL genome for given release/species.

    Parameters
    ----------
    species : str
        Species latin name.
    release : int
        Release number. Only ENSEMBL releases {0}-{1} are available.
    out_dir : str
        Download to this directory (if not given, current working directory).
    genome : str
        Genome filename (must have .gz file extension). If not given,
        species.release.fa.gz is used. If genome is provided as absolute path,
        value of out_dir parameter is ignored and file is saved to given
        absolute path.
    chromosomes : list_str
        If given, do not download the whole genome, but listed
        chromosomes only. Chromosomes can be given as strings or integers.

    Returns
    -------
    str
        Downloaded genome/sequnce filename.

    """
    iCount.log_inputs(LOGGER, level=logging.INFO)

    if species not in iCount.genomes.ensembl.species():
        raise ValueError('Invalid species name.')

    if not MIN_RELEASE_SUPPORTED <= release <= MAX_RELEASE_SUPPORTED:
        raise ValueError('Release should be a number between {} and {}.'.format(
            MIN_RELEASE_SUPPORTED, MAX_RELEASE_SUPPORTED))

    if genome:
        assert genome.endswith(('.fa', '.fasta', '.fa.gz', '.fasta.gz'))

    # If absolute path is given, ignore out_dir:
    if genome and os.path.isabs(genome):
        if out_dir:
            LOGGER.info('out_dir parameter has been changed, since absolute '
                        'path is provided by genome parameter.')
        out_dir = os.path.dirname(genome)
        genome = os.path.basename(genome)
    if not out_dir:
        out_dir = os.getcwd()
    if not os.path.isdir(out_dir):
        raise ValueError('Directory "{}" does not exist.'.format(out_dir))

    ftp = iCount.genomes.get_ftp_instance(BASE_URL)
    server_dir = '/pub/release-{}/fasta/{}/dna'.format(release, species)
    ftp.cwd(server_dir)
    all_fasta_files = ftp.nlst()

    filtered_files = []

    # Recognize all "whole-chromosme" files
    regex = r'{}\.[\w\d\-]+\.(\d+\.)*dna\.chromosome\.' \
            r'[a-zA-Z0-9]+\.fa\.gz'.format(species.capitalize())
    for fasta_file in all_fasta_files:
        if re.match(regex, fasta_file):
            filtered_files.append(fasta_file)

    def sorting_func(name):
        """Sort names by chromosome order."""
        key = name.split('.')[-3]
        if key == 'MT':
            return 'ZZ'  # this makes mitohondrial "chromosome" the last one
        return key
    # Sort the filtered_files in correct order of chromosomes:
    filtered_files.sort(key=sorting_func)

    # If parameter chromosomes is given, take only a subset of all chromosomes
    if chromosomes:
        subset_list = []
        for chromosome in chromosomes:
            regex = r'.*chromosome\.{}\.fa\.gz'.format(str(chromosome))
            chromosome_found = False
            for file_ in filtered_files:
                if re.match(regex, file_):
                    subset_list.append(file_)
                    chromosome_found = True
            if not chromosome_found:
                # Provide user with a list of available chromosomes:
                chrom_regex = r'.*chromosome\.([a-zA-Z0-9]+)\.fa\.gz'
                available = [re.match(chrom_regex, file_).group(1) for file_ in filtered_files]
                raise ValueError('Could not find chromosome {}. Available chromosomes '
                                 'are: {}'.format(chromosome, ' '.join(available)))
        filtered_files = subset_list
        if not genome:
            genome = '{}.{}.chr{}.fa.gz'.format(
                species, release, '_'.join(map(str, chromosomes)))

    if not genome:
        genome = '{}.{}.fa.gz'.format(species, release)
    target_path = os.path.abspath(os.path.join(out_dir, genome))

    LOGGER.info('Downloading FASTA file into: %s', target_path)
    temp_dir = os.path.join(iCount.TMP_ROOT, 'ensembl')
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    for fname in filtered_files:
        LOGGER.debug('Downloading file: %s', fname)
        # Download all files to tempdir:
        temp_file = os.path.join(temp_dir, fname)
        with open(temp_file, 'wb') as fhandle:
            ftp.retrbinary('RETR ' + fname, fhandle.write)

        # Write the content into final file (target_fname)
        with gzip.open(temp_file, 'rb') as f_in, gzip.open(target_path, 'ab') as f_out:
            shutil.copyfileobj(f_in, f_out)

    # Clean up: delete temporary files:
    ftp.quit()
    shutil.rmtree(temp_dir)

    # Compute chromosome lengths:
    chrom_length(target_path)

    LOGGER.info('Done.')
    return target_path
