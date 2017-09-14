""".. Line to protect from pydocstyle D205, D400.

GENCODE API
-----------

Functions to query and download data from the `Gencode`_ FTP site.

.. _Gencode:
    https://www.gencodegenes.org/

"""

import logging
import os
import re

import iCount
from iCount.genomes.ensembl import get_ftp_instance

BASE_URL = 'ftp.sanger.ac.uk'
LOGGER = logging.getLogger(__name__)


def _to_int(string):
    """Convert string to integer, if possible."""
    try:
        return int(string)
    except ValueError:
        return string


def species():
    """
    Get list of available species.

    Returns
    -------
    list
        List of species.

    """
    iCount.log_inputs(LOGGER, level=logging.INFO)

    ftp = get_ftp_instance(BASE_URL)
    ftp.cwd('pub/gencode')
    spec_list = [item[8:] for item in ftp.nlst() if re.match(r'Gencode_\w+', item)]
    ftp.quit()
    LOGGER.info('There are %d species available: %s',
                len(spec_list), ','.join(map(str, spec_list)))
    return spec_list


# pylint: disable=redefined-outer-name
def releases(species='human'):
    """
    Get list of available releases for given species.

    Parameters
    ----------
    species : str
        Species name.

    Returns
    -------
    list
        List of available releases

    """
    if species is None:
        species = 'human'
    ftp = get_ftp_instance(BASE_URL)
    # set current working directory
    ftp.cwd('pub/gencode/Gencode_{}'.format(species))

    fetch = [item[8:] for item in ftp.nlst() if re.match(r'release_\w+', item)]
    ftp.quit()

    def _keyfunc(item):
        """Separate mixed alphanumeric string into numbers and words."""
        return [_to_int(i) for i in re.match(r'([0-9]*)([a-zA-Z]*)', item).groups()]
    out = sorted(fetch, reverse=True, key=_keyfunc)

    LOGGER.info('There are %d GENCODE releases available for %s: %s',
                len(out), species, ','.join(map(str, (out))))
    return out


# pylint: disable=redefined-outer-name
def annotation(species, release, out_dir=None, annotation=None):
    """
    Download GENCODE annotation for given release/species.

    Parameters
    ----------
    species : str
        Species name.
    release : int
        Release number.
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

    if species not in iCount.genomes.gencode.species():
        raise ValueError('Invalid species name.')

    if str(release) not in iCount.genomes.gencode.releases(species):
        raise ValueError('Invalid release number.')

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
        raise ValueError('Directory "{}" does not exist'.format(out_dir))

    ftp = get_ftp_instance(BASE_URL)
    ftp.cwd('/pub/gencode/Gencode_{}/release_{}/'.format(species, release))
    server_files = ftp.nlst()

    # Get the desired annotation file form list of all server files:
    regex = r'gencode\.v{}\.annotation\.gtf\.gz'.format(release)
    annotation_file = next((fname for fname in server_files if re.match(regex, fname)), None)

    if not annotation_file:
        LOGGER.info('No GTF file found for species %s, release %s', species, str(release))
        return None

    saveas_fname = os.path.join(out_dir, annotation)
    LOGGER.info('Downloading GTF to: %s', saveas_fname)
    with open(saveas_fname, 'wb') as fhandle:
        ftp.retrbinary('RETR ' + annotation_file, fhandle.write)
    ftp.quit()

    LOGGER.info('Done.')
    return os.path.abspath(saveas_fname)


# pylint: disable=redefined-outer-name
def genome(species, release, out_dir=None, genome=None):
    """
    Download GENCODE genome for given release/species.

    Parameters
    ----------
    species : str by using samtools faidx
        Species latin name.
    release : int
        Release number. Only gencode releases {0}-{1} are available.
    out_dir : str
        Download to this directory (if not given, current working directory).
    genome : str
        Genome filename (must have .gz file extension). If not given,
        species.release.fa.gz is used. If genome is provided as absolute path,
        value of out_dir parameter is ignored and file is saved to given
        absolute path.

    Returns
    -------
    str
        Downloaded genome/sequnce filename.

    """
    iCount.log_inputs(LOGGER, level=logging.INFO)

    if species not in iCount.genomes.gencode.species():
        raise ValueError('Invalid species name.')

    if str(release) not in iCount.genomes.gencode.releases(species):
        raise ValueError('Invalid release number.')

    if genome:
        assert genome.endswith(('.fa', '.fasta', '.fa.gz', '.fasta.gz'))
    else:
        genome = '{}.{}.fa.gz'.format(species, release)

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
        raise ValueError('Directory "{}" does not exist'.format(out_dir))

    ftp = get_ftp_instance(BASE_URL)
    ftp.cwd('/pub/gencode/Gencode_{}/release_{}/'.format(species, release))
    regex = r'[\w\d\.]+\.genome\.fa\.gz'
    fasta_file = next((fname for fname in ftp.nlst() if
                       (re.match(regex, fname) and 'primary_assembly' not in fname)), None)

    target_path = os.path.abspath(os.path.join(out_dir, genome))
    LOGGER.info('Downloading FASTA file into: %s', target_path)
    with open(target_path, 'wb') as fhandle:
        ftp.retrbinary('RETR ' + fasta_file, fhandle.write)
    ftp.quit()

    # Compute chromosome lengths:
    iCount.genomes.ensembl.chrom_length(target_path)

    LOGGER.info('Done.')
    return target_path
