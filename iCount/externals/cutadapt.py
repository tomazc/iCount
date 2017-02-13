""".. Line to protect from pydocstyle D205, D400.

Cutadapt
--------

Remove adapter sequences from reads in FASTQ file.
"""

import subprocess


def get_version():
    """Get cutadapt version."""
    args = ['cutadapt', '--version']
    try:
        ver = subprocess.check_output(args, shell=False, universal_newlines=True)
        return str(ver).rstrip('\n\r')
    except (FileNotFoundError, subprocess.CalledProcessError):
        return None


def run(reads, reads_trimmed, adapter, qual_base=64, qual_trim=None, minimum_length=None):
    """
    Remove adapter sequences from high-throughput sequencing reads.

    Parameters
    ----------
    reads : str
        Input FASTQ file.
    reads_trimmed : str
        Output FASTQ file containing trimmed reads.
    adapter : str
        Sequence of an adapter ligated to the 3' end.
    qual_base : int
        Assume that quality values in FASTQ are encoded as ascii(quality +
        QUALITY_BASE). value of 64 corresponds to old Illumina FASTQ files.
    qual_trim : int
        Trim low-quality bases before adapter removal.
    minimum_length : int
        Discard trimmed reads that are shorter than `minimum_length`.

    Returns
    -------
    int
        Return code of the `cutadapt` program.

    """
    args = [
        'cutadapt', '--quiet',
        '-a', adapter,
        '--quality-base={:d}'.format(qual_base),
    ]
    if qual_trim is not None:
        args.extend(['-q', '{:d}'.format(qual_trim)])
    if minimum_length is not None:
        args.extend(['-m', '{:d}'.format(minimum_length)])
    args.extend(['-o', reads_trimmed, reads])

    return subprocess.call(args, shell=False)
