""".. Line to protect from pydocstyle D205, D400.

Files
=====

iCount works with various formats that store `FASTA`_ and `FASTQ`_ sequencing data, `GTF`_ genome
annotation, `BAM`_ data on mapped reads, `BED`_ files with quantified cross-linked sites.

.. autofunction:: iCount.files.gz_open
.. autofunction:: iCount.files.decompress_to_tempfile

.. automodule:: iCount.files.bed
   :members:

.. automodule:: iCount.files.fastq
   :members:

.. automodule:: iCount.files.fasta
   :members:

.. automodule:: iCount.files.gtf
   :members:

.. _FASTA:
    https://en.wikipedia.org/wiki/FASTA_format

.. _FASTQ:
    https://en.wikipedia.org/wiki/FASTQ_format

.. _GTF:
    http://www.gencodegenes.org/gencodeformat.html

.. _BAM:
    https://samtools.github.io/hts-specs/SAMv1.pdf

.. _BED:
    http://bedtools.readthedocs.io/en/latest/content/general-usage.html#bed-format

"""

import os
import gzip
import tempfile
import shutil

import iCount

from . import bed
from . import fasta
from . import fastq


def gz_open(fname, omode):
    """
    Use :py:mod:`gzip` library to open compressed files ending with .gz.

    Parameters
    ----------
    fname : str
        Path to file to open.
    omode : str
        String indicating how the file is to be opened.

    Returns
    -------
    file
        File Object.

    """
    if fname.endswith(".gz"):
        return gzip.open(fname, omode)
    return open(fname, omode)


def decompress_to_tempfile(fname, context='misc'):
    """
    Decompress files ending with .gz to a temporary file and return filename.

    If file does nto end with .gz, juts return fname.

    Parameters
    ----------
    fname : str
        Path to file to open.
    context : str
        Name of temporary subfolder where temporary file is created.

    Returns
    -------
    str
        Path to decompressed file.

    """
    if fname.endswith('.gz'):
        tmp_dir = os.path.join(iCount.TMP_ROOT, context)
        if not os.path.exists(tmp_dir):
            os.makedirs(tmp_dir)

        suffix = '_{:s}'.format(os.path.basename(fname))
        fout = tempfile.NamedTemporaryFile(suffix=suffix, dir=tmp_dir, delete=False)
        fin = gzip.open(fname, 'r')
        shutil.copyfileobj(fin, fout)
        fin.close()
        fout.close()
        return fout.name

    return fname


def _f2s(number, dec=4):
    """
    Return string representation of ``number``.

    Returned string is:

        * without trailing decimal zeros,
        * with at most ``dec`` decimal places.
    """
    if not isinstance(number, (int, float)):
        return number
    return '{{:.{:d}f}}'.format(dec).format(number).rstrip('0').rstrip('.')
