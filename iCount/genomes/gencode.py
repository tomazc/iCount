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

BASE_URL = 'ftp.sanger.ac.uk'
LOGGER = logging.getLogger(__name__)


MIN_RELEASE_SUPPORTED = {
    'human': '22',
    'mouse': 'M5',
}


def species():
    """
    Get list of available species.

    Returns
    -------
    list
        List of species.

    """
    iCount.log_inputs(LOGGER, level=logging.INFO)

    ftp = iCount.genomes.get_ftp_instance(BASE_URL)
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
    ftp = iCount.genomes.get_ftp_instance(BASE_URL)
    # set current working directory
    ftp.cwd('pub/gencode/Gencode_{}'.format(species))

    fetch = [item[8:] for item in ftp.nlst() if re.match(r'release_\w+', item)]
    ftp.quit()

    def _keyfunc(item):
        """Separate mixed alphanumeric string into numbers and words."""
        # pylint: disable=protected-access
        regex = r'([a-zA-Z]*)([0-9]*)([a-zA-Z]*)'
        return [iCount.genomes._to_int(i) for i in re.match(regex, item).groups()]
    out = sorted(fetch, reverse=True, key=_keyfunc)
    # Suggest only supported releases: the ones with higher version than in MIN_RELEASE_SUPPORTED
    until_here = out.index(MIN_RELEASE_SUPPORTED[species])
    out = out[:until_here + 1]

    LOGGER.info('There are %d GENCODE releases available for %s: %s',
                len(out), species, ','.join(map(str, (out))))
    return out


# pylint: disable=redefined-outer-name
def annotation(species, release, out_dir=None, annotation=None):
    """
    Download GENCODE annotation for given release/species.

    Note: This will download the "primary assembly" type of annotation.

    Parameters
    ----------
    species : str
        Species name.
    release : str
        Release number.
    out_dir : str
        Download to this directory (if not given, current working directory).
    annotation : str
        Annotation filename (must have .gz file extension). If not given
        original filename will be used. If annotation is provided as absolute
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

    # If absolute path is given, ignore out_dir:
    if annotation is not None and os.path.isabs(annotation):
        if out_dir:
            LOGGER.info('out_dir parameter has been changed, since absolute '
                        'path is provided by annotation parameter.')
        out_dir = os.path.dirname(annotation)
        annotation = os.path.basename(annotation)
    if not out_dir:
        out_dir = os.getcwd()
    if not os.path.isdir(out_dir):
        raise ValueError('Directory "{}" does not exist'.format(out_dir))

    ftp = iCount.genomes.get_ftp_instance(BASE_URL)
    ftp.cwd('/pub/gencode/Gencode_{}/release_{}/'.format(species, release))
    server_files = ftp.nlst()

    # Get the desired annotation file form list of all server files:
    regex = r'gencode\.v{}\.primary_assembly\.annotation\.gtf\.gz'.format(release)
    annotation_file = next((fname for fname in server_files if re.match(regex, fname)), None)

    if not annotation_file:
        LOGGER.info('No GTF file found for species %s, release %s', species, str(release))
        return None

    target_path = os.path.join(out_dir, annotation or annotation_file)
    LOGGER.info('Downloading GTF to: %s', target_path)
    with open(target_path, 'wb') as fhandle:
        ftp.retrbinary('RETR ' + annotation_file, fhandle.write)
    ftp.quit()

    LOGGER.info('Done.')
    return os.path.abspath(target_path)


# pylint: disable=redefined-outer-name
def genome(species, release, out_dir=None, genome=None):
    """
    Download GENCODE genome for given release/species.

    Note: This will download the "primary assembly" type of genome.

    Parameters
    ----------
    species : str by using samtools faidx
        Species latin name.
    release : str
        Release number. Only gencode releases {0}-{1} are available.
    out_dir : str
        Download to this directory (if not given, current working directory).
    genome : str
        Genome filename (must have .gz file extension). If not given original
        filename will be used. If genome is provided as absolute path, value
        of out_dir parameter is ignored and file is saved to given absolute
        path.

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

    # If absolute path is given, ignore out_dir:
    if genome is not None and os.path.isabs(genome):
        if out_dir:
            LOGGER.info('out_dir parameter has been changed, since absolute '
                        'path is provided by genome parameter.')
        out_dir = os.path.dirname(genome)
        genome = os.path.basename(genome)
    if not out_dir:
        out_dir = os.getcwd()
    if not os.path.isdir(out_dir):
        raise ValueError('Directory "{}" does not exist'.format(out_dir))

    ftp = iCount.genomes.get_ftp_instance(BASE_URL)
    ftp.cwd('/pub/gencode/Gencode_{}/release_{}/'.format(species, release))
    regex = r'[\w\d\.]+\.primary_assembly\.genome\.fa\.gz'
    fasta_file = next((fname for fname in ftp.nlst() if re.match(regex, fname)), None)

    target_path = os.path.abspath(os.path.join(out_dir, genome or fasta_file))
    LOGGER.info('Downloading FASTA file into: %s', target_path)
    with open(target_path, 'wb') as fhandle:
        ftp.retrbinary('RETR ' + fasta_file, fhandle.write)
    ftp.quit()

    # Compute chromosome lengths:
    iCount.genomes.ensembl.chrom_length(target_path)

    LOGGER.info('Done.')
    return target_path
