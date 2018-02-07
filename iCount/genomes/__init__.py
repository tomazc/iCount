""".. Line to protect from pydocstyle D205, D400.

Genomes
=======

This module provides access to `Ensembl`_ and `Gencode`_ genome sequence and annotation.
Also, segmentation into genes and segments of same type (exon, intron, UTR, ...) is supported.

.. automodule:: iCount.genomes.ensembl
   :members:

.. automodule:: iCount.genomes.gencode
   :members:

.. automodule:: iCount.genomes.segment
   :members:

.. _Ensembl:
    http://www.ensembl.org/index.html

.. _Gencode:
    https://www.gencodegenes.org/

"""
import ftplib

from . import ensembl
from . import gencode
from . import segment


SUPPORTED_SOURCES = [
    'ensembl',
    'gencode',
]


def _to_int(string):
    """Convert string to integer, if possible."""
    try:
        return int(string)
    except (ValueError, TypeError):
        return string


def get_ftp_instance(base_url):
    """
    Get ftplib.FTP object that is connected to base_url.

    Returns
    -------
    ftplib.FTP
        FTP object connected to base_url.

    """
    try:
        ftp = ftplib.FTP(base_url)
        ftp.login()
        return ftp
    except Exception:
        raise ConnectionError('Problems connecting to ENSEMBL FTP server.')


# pylint: disable=redefined-outer-name
def species(source='gencode', release=None):
    """
    Get list of available species.

    Parameters
    ----------
    source : str
        Source of data. Only ENSEMBL or GENCODE are available.
    release : str
        Release number. Only relevant if source is ENSEMBL.

    Returns
    -------
    list
        List of species.

    """
    source = source.lower()
    if source not in SUPPORTED_SOURCES:
        raise ValueError('Source {} is not supported.'.format(source))

    if source == 'gencode':
        return gencode.species()
    else:
        return ensembl.species(release=_to_int(release))  # pylint:disable=protected-access


def releases(source='gencode', species=None):
    """
    Get list of available releases.

    Parameters
    ----------
    source : str
        Source of data. Only ENSEMBL or GENCODE are available.
    species : str
        Species name. Only relevant if source is GENCODE.

    Returns
    -------
    list
        List of available releases

    """
    source = source.lower()
    if source not in SUPPORTED_SOURCES:
        raise ValueError('Source {} is not supported.'.format(source))

    if source == 'gencode':
        return gencode.releases(species=species)
    else:
        return ensembl.releases()


# pylint: disable=redefined-outer-name,protected-access
def annotation(species, release, out_dir=None, annotation=None, source='gencode'):
    """
    Download annotation for given release/species/source.

    Parameters
    ----------
    species : str
        Species name.
    release : str
        Release number.
    out_dir : str
        Download to this directory (if not given, current working directory).
    annotation : str
        Annotation filename (must have .gz file extension). If annotation is
        provided as absolute path, value of out_dir parameter is ignored and
        file is saved to given absolute path.
    source : str
        Source of data. Only ENSEMBL or GENCODE are available.

    Returns
    -------
    str
        Downloaded annotation filename.

    """
    source = source.lower()
    if source not in SUPPORTED_SOURCES:
        raise ValueError('Source {} is not supported.'.format(source))

    if source == 'gencode':
        return gencode.annotation(
            species=species, release=release, out_dir=out_dir, annotation=annotation)
    else:
        return ensembl.annotation(species=species, release=_to_int(release),
                                  out_dir=out_dir, annotation=annotation)


# pylint: disable=redefined-outer-name,protected-access
def genome(species, release, out_dir=None, genome=None, chromosomes=None, source='gencode'):
    """
    Download genome for given release/species/source.

    Parameters
    ----------
    species : str
        Species name.
    release : str
        Release number.
    out_dir : str
        Download to this directory (if not given, current working directory).
    genome : str
        Genome filename (must have .gz file extension). If genome is provided
        as absolute path, value of out_dir parameter is ignored and file is
        saved to given absolute path.
    chromosomes : list_str
        If given, do not download the whole genome, but listed
        chromosomes only. Only relevant if source is ENSEMBL.
    source : str
        Source of data. Only ENSEMBL or GENCODE are available.

    Returns
    -------
    str
        Downloaded genome/sequnce filename.

    """
    source = source.lower()
    if source not in SUPPORTED_SOURCES:
        raise ValueError('Source {} is not supported.'.format(source))

    if source == 'gencode':
        return gencode.genome(
            species=species, release=release, out_dir=out_dir, genome=genome)
    else:
        return ensembl.genome(species=species, release=_to_int(release), out_dir=out_dir,
                              genome=genome, chromosomes=chromosomes)
