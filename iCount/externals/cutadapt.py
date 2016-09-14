"""
Interface to running cutadapt.
"""

import subprocess


def get_version():
    args = [
        'cutadapt', '--version',
    ]
    try:
        ver = subprocess.check_output(args, shell=False,
                                      universal_newlines=True)
        return str(ver).rstrip('\n\r')
    except:
        return None


def run(in_fname, out_fname, adapter, qual_base=64, qual_trim=None,
        minimum_length=None):
    args = [
        'cutadapt', '--quiet',
        '-a', adapter,
        '--quality-base={:d}'.format(qual_base),
    ]
    if qual_trim is not None:
        args.extend(['-q', '{:d}'.format(qual_trim)])
    if minimum_length is not None:
        args.extend(['-m', '{:d}'.format(minimum_length)])
    args.extend(['-o', out_fname, in_fname])

    return subprocess.call(args, shell=False)
