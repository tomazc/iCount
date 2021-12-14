""".. Line to protect from pydocstyle D205, D400.

Cutadapt
--------

Remove adapter sequences from reads in FASTQ file.
"""

import os
import shutil
import subprocess
import tempfile

import iCount
from iCount.files.fastq import get_qual_encoding, ENCODING_TO_OFFSET


def get_version():
    """Get cutadapt version."""
    args = ['cutadapt', '--version']
    try:
        ver = subprocess.check_output(args, shell=False, universal_newlines=True)
        return str(ver).rstrip('\n\r')
    except (FileNotFoundError, subprocess.CalledProcessError):
        return None


def convert_version(version, n=3):
    """Converts string representation of a version into a comparable sematic 
    version format. N can be used to specify how far it should check. By default,
    N is set to 3. This corresponds to the sematic versioning specifications of
    'MAJOR.MINOR.PATCH'.

    Parameters
    ----------
    version : str
        Cutadapt version that is in the user's PATH. Example: 1.14

    Returns
    -------
    version_mask: tuple(int, int, int)
        Tuple containing integers corresponding to MAJOR, MINOR, PATCH 
        components of a given semantic version.
    """
    # Tuples containing integers can
    # be directly compared in python. 
    # As an example:
    # (1, 2, 0) > (1, 1, 0)
    # returns True 
    version_mask = [0] * n # Pad array with N zeros
    if version is None:
        return tuple(version_mask)

    for i,v in enumerate(version.split('.')):
        if v.startswith('v'):
            # Remove any characters starting with v
            # Example: v1 -> 1
            v = v.replace('v', '')
        if i < n:
            # Check up to 3 semantic version components,
            # i.e. MAJOR.MINOR.PATCH
            version_mask[i] = int(v)

    version_mask = tuple(version_mask)

    return version_mask


def multithreading_supported(version, min_version='1.15'):
    """Checks the version of cutadapt to see if multithreading is supported.
    Older versions of cutadapt do NOT support multi-treading. This feature was
    added in cutadapt version '1.15'. If the version of cutadapt supports 
    multithreading and the user provides a value to --threads option of the 
    demultiplex sub command, then cutadapt will be run with the -j option.

    Parameters
    ----------
    version : str
        Cutadapt version that is in the user's PATH.
    min_version : str
        Minimum version that of cutadapt supports multithreading 

    Returns
    -------
    boolean
        Return whether cutadapt supports multithreading, where
        True indicates that multi-threading is support. 
    """
    # Default to not supporting multithreading 
    supported = False

    # Check if user version is greater than 
    # or equal to version '1.15.0'. If so,
    # multi-threading is supported. 
    min_sematic_version = convert_version(min_version)
    user_sematic_version = convert_version(version)
    if user_sematic_version:
        supported = user_sematic_version >= min_sematic_version

    return supported


def run(reads, adapter, reads_trimmed=None, overwrite=False, qual_trim=None, minimum_length=None, overlap=None,
        untrimmed_output=None, error_rate=None, threads=1):
    """
    Remove adapter sequences from high-throughput sequencing reads.

    Parameters
    ----------
    reads : str
        Input FASTQ file.
    adapter : str
        Sequence of an adapter ligated to the 3' end.
    reads_trimmed : str
        Output FASTQ file containing trimmed reads. If not provided
    overwrite : bool
        If true, overwrite input file (reads) with trimmed file.
    qual_trim : int
        Trim low-quality bases before adapter removal.
    minimum_length : int
        Discard trimmed reads that are shorter than `minimum_length`.
    overlap : int
        Require ``overlap`` overlap between read and adapter for an
        adapter to be found.
    untrimmed_output : str
        Write reads that do not contain any adapter to this file.
    error_rate : float
        Maximum allowed error rate (no. of errors divided by the length
        of the matching region).
    threads : int
        Number of CPU cores to use with cutadapt. This feature is only enabled with
        versions of cutadapt greater than or equal to 1.15.

    Returns
    -------
    int
        Return code of the `cutadapt` program.

    """
    args = [
        'cutadapt',
        '--quiet',
        '-a', adapter,
    ]
    qual_base = ENCODING_TO_OFFSET.get(get_qual_encoding(reads), 33)
    args.extend(['--quality-base={}'.format(qual_base)])

    if reads_trimmed is None:
        # Auto-generate output name:
        extension = '.gz' if reads.endswith('.gz') else ''
        name = next(tempfile._get_candidate_names()) + '.fq' + extension  # pylint: disable=protected-access
        reads_trimmed = os.path.join(iCount.TMP_ROOT, name)
    if qual_trim is not None:
        args.extend(['-q', '{:d}'.format(qual_trim)])
    if minimum_length is not None:
        args.extend(['-m', '{:d}'.format(minimum_length)])
    if overlap is not None:
        args.extend(['--overlap', '{:d}'.format(overlap)])
    if untrimmed_output is not None:
        args.extend(['--untrimmed-output', '{}'.format(untrimmed_output)])
    if error_rate is not None:
        args.extend(['--error-rate', '{}'.format(error_rate)])
    if multithreading_supported(get_version()):
        args.extend(['-j', '{}'.format(threads)])
    args.extend(['-o', reads_trimmed, reads])

    rcode = subprocess.call(args, shell=False)

    if overwrite:
        shutil.move(reads_trimmed, reads)

    return rcode


if __name__ == '__main__':
    """Unit-testing"""
    # Testing functionality of convert_version() 
    assert convert_version('v1.15', 2) == (1, 15)
    assert convert_version('1.15', 3) == (1, 15, 0)
    assert convert_version(None) == (0, 0, 0)

    # Testing functionality of multithreading_supported()
    assert multithreading_supported('1.15') == True
    assert multithreading_supported('1.14') == False 
    assert multithreading_supported(None) == False 

    # Mocking integration of convert_version(), 
    # multithreading_supported(), get_version(),
    # within run()
    if multithreading_supported(get_version()):
        print("Adding... -j option to cutadapt command")