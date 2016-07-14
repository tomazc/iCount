"""
Interface to running cutadapt.
"""

import subprocess


def get_version():
    args = ['cutadapt', '--version']
    try:
        ver = subprocess.check_output(args, shell=False,
                                      universal_newlines=True)
        return str(ver).rstrip('\n\r')
    except:
        return None


def run(in_fname, out_fname, adapter, qual_base=64, qual_trim=None):
    args = ['cutadapt', '--quiet', '-a', adapter,
            '--quality-base=%s' % qual_base]
    if qual_trim is not None:
        args.extend(['-q', '%s' % qual_trim])
    args.extend(['-o', out_fname, in_fname])

    return subprocess.call(args, shell=False)
