"""
Interface to running STAR.
"""

import os
import subprocess

import iCount


def get_version():
    args = ['STAR', '--version']
    try:
        ver = subprocess.check_output(args, shell=False,
                                      universal_newlines=True)
        return str(ver).rstrip('\n\r')
    except:
        return None


def build_index(genome_fname, outdir, annotation='', overhang=100, threads=1):
    genome_fname2 = iCount.files.temp.decompress_to_tempfile(
        genome_fname, 'starindex')
    args = [
        'STAR',
        '--runThreadN', '{:d}'.format(threads),
        '--runMode', 'genomeGenerate',
        '--genomeDir', '{:s}'.format(outdir),
        '--genomeFastaFiles', '{:s}'.format(genome_fname2),
        '--sjdbOverhang', '{:d}'.format(overhang),
    ]
    if annotation:
        annotation2 = iCount.files.temp.decompress_to_tempfile(annotation,
                                                               'starindex')
        args.extend([
            '--sjdbGTFfile', annotation2,
        ])
    else:
        annotation2 = annotation

    try:
        ret_code = subprocess.call(args, shell=False)
    finally:
        # remove temporary decompressed files
        if genome_fname != genome_fname2:
            os.remove(genome_fname2)
        if annotation != annotation2:
            os.remove(annotation2)

    return ret_code


def map_reads(sequences_fname, genomedir, outdir, annotation='',
              multimax=10, mismatches=2, threads=1):
    sequences_fname2 = iCount.files.temp.decompress_to_tempfile(
        sequences_fname, 'starmap')

    args = [
        'STAR',
        '--runMode', 'alignReads',
        '--runThreadN', '{:d}'.format(threads),
        '--genomeDir', '{:s}'.format(genomedir),
        '--readFilesIn', '{:s}'.format(sequences_fname2),
    ]
    if not outdir.endswith('/'):
        outdir += '/'

    args.extend([
        '--outFileNamePrefix', '{:s}'.format(outdir),
        '--outSAMprimaryFlag', 'AllBestScore',
        '--outFilterMultimapNmax', '{:d}'.format(multimax),
        '--outFilterMismatchNmax', '{:d}'.format(mismatches),
        '--alignEndsType', 'EndToEnd',
        # otherwise soft-clipping of the starts and ends may produce too
        # many multiple hits
    ])
    if annotation:
        annotation2 = iCount.files.temp.decompress_to_tempfile(annotation,
                                                               'starmap')
        args.extend([
            '--sjdbGTFfile', annotation,
        ])
    else:
        annotation2 = annotation

    try:
        ret_code = subprocess.call(args, shell=False)
    finally:
        # remove temporary decompressed files
        if sequences_fname != sequences_fname2:
            os.remove(sequences_fname2)
        if annotation != annotation2:
            os.remove(annotation2)

    return ret_code
