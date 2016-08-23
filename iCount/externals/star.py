"""
Interface to running STAR.
"""

import subprocess


def get_version():
    args = ['STAR', '--version']
    try:
        ver = subprocess.check_output(args, shell=False,
                                      universal_newlines=True)
        return str(ver).rstrip('\n\r')
    except:
        return None


def run(genome_fname, outdir, annotation='', overhang=100, threads=1):
    args = ['STAR',
            '--runThreadN', '{:d}'.format(threads),
            '--runMode', 'genomeGenerate',
            '--genomeDir', '{:s}'.format(outdir),
            '--genomeFastaFiles', '{:s}'.format(genome_fname),
            '--sjdbOverhang', '{:d}'.format(overhang),
            ]
    if annotation:
        args.extend([
            '--sjdbGTFfile', annotation,
        ])

    return subprocess.call(args, shell=False)
