""".. Line to protect from pydocstyle D205, D400.

STAR aligner
------------

Interface to running STAR.
"""

import os
import logging
import subprocess

from threading import Thread
from queue import Queue

import iCount

LOGGER = logging.getLogger(__name__)


def _execute(cmd):
    """Give stdout and stderr in realtime when executing the cmd as iterator."""
    def readline_output(out, queue, name):  # pylint: disable=missing-docstring
        for line in iter(out.readline, ''):
            queue.put((name, line))
        out.close()
        queue.put((name, 'readline_output finished.'))

    popen = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             universal_newlines=True)

    queue = Queue()
    Thread(target=readline_output, args=(popen.stdout, queue, 'stdout'), daemon=True).start()
    Thread(target=readline_output, args=(popen.stderr, queue, 'stderr'), daemon=True).start()

    done = 0
    while True:
        out, message = queue.get()
        if message == 'readline_output finished.':
            done += 1
        else:
            yield ('{}_line'.format(out), message)

        if done >= 2:
            break

    yield ('return_code', popen.wait())


def get_version():
    """Get STAR version."""
    args = ['STAR', '--version']
    try:
        ver = subprocess.check_output(args, shell=False,
                                      universal_newlines=True)
        return str(ver).rstrip('\n\r')
    except (FileNotFoundError, subprocess.CalledProcessError):
        return None


def build_index(genome, genome_index, annotation='', overhang=100, overhang_min=8, threads=1,
                genome_sasparsed=1, genome_saindexnbases=14):
    """
    Call STAR to generate genome index, which is used for mapping.

    Parameters
    ----------
    genome : str
        Genome sequence to index.
    genome_index : str
        Output folder, where to store genome index.
    annotation : str
        Annotation that defines splice junctions.
    overhang : int
        Sequence length around annotated junctions to be used by STAR when
        constructing splice junction database.
    overhang_min : int
        Minimum overhang for unannotated junctions.
    threads : int
        Number of threads that STAR can use for generating index.
    genome_sasparsed : int
        STAR parameter genomeSAsparseD.
        Suffix array sparsity. Bigger numbers decrease RAM requirements
        at the cost of mapping speed reduction. Suggested values
        are 1 (30 GB RAM) or 2 (16 GB RAM).
    genome_saindexnbases : int
        STAR parameter genomeSAindexNbases.
        SA pre-indexing string length, typically between 10 and 15.
        Longer strings require more memory, but result in faster searches.

    Returns
    -------
    int
        Star return code.

    """
    iCount.log_inputs(LOGGER, level=logging.INFO)

    if not os.path.isdir(genome_index):
        raise FileNotFoundError('Output directory does not exist. Make sure it does.')

    LOGGER.info('Building genome index with STAR for genome %s', genome)
    genome_fname2 = iCount.files.decompress_to_tempfile(
        genome, 'starindex')
    args = [
        'STAR',
        '--runThreadN', '{:d}'.format(threads),
        '--genomeSAsparseD', '{:d}'.format(genome_sasparsed),
        '--genomeSAindexNbases', '{:d}'.format(genome_saindexnbases),
        '--runMode', 'genomeGenerate',
        '--genomeDir', '{:s}'.format(genome_index),
        '--genomeFastaFiles', '{:s}'.format(genome_fname2),
        '--alignSJoverhangMin', '{:d}'.format(overhang_min),
        '--outFileNamePrefix', '{:s}'.format(genome_index),
    ]
    if annotation:
        annotation2 = iCount.files.decompress_to_tempfile(annotation, 'starindex')
        args.extend([
            '--sjdbGTFfile', annotation2,
            '--sjdbOverhang', '{:d}'.format(overhang),
        ])
    else:
        annotation2 = annotation

    try:
        ret_code = 1
        for name, value in _execute(args):
            if name == 'return_code' and isinstance(value, int):
                ret_code = value
                break
            elif name == 'stdout_line':
                LOGGER.info(value.strip())
            elif name == 'stderr_lines':
                for line in value.split('\n'):
                    LOGGER.error(line.strip())
    finally:
        # remove temporary decompressed files
        if genome != genome_fname2:
            os.remove(genome_fname2)
        if annotation != annotation2:
            os.remove(annotation2)

    LOGGER.info('Done.')
    return ret_code


def map_reads(reads, genome_index, out_dir, annotation='', multimax=10, mismatches=2, threads=1, genome_load=False):
    """
    Map FASTQ file reads to reference genome.

    Parameters
    ----------
    reads : str
        Sequencing reads to map to genome.
    genome_index : str
        Folder with genome index.
    out_dir : str
        Output folder, where to store mapping results.
    annotation : str
        GTF annotation needed for mapping to splice-junctions.
    multimax : int
        Number of allowed multiple hits.
    mismatches : int
        Number of allowed mismatches.
    threads : int
        Number of threads that STAR can use for generating index.
    genome_load: bool
        Load genome into shared memory.
        Shared memory must be available in the system. See Chapter 3.3 in STAR manual.

    Returns
    -------
    int
        Return code

    """
    iCount.log_inputs(LOGGER, level=logging.INFO)

    if not os.path.isdir(genome_index):
        raise FileNotFoundError('Directory with genome index does not exist. Make sure it does.')
    if not os.path.isdir(out_dir):
        raise FileNotFoundError('Output directory does not exist. Make sure it does.')

    LOGGER.info('Mapping reads from %s', reads)
    sequences_fname2 = iCount.files.decompress_to_tempfile(
        reads, 'starmap')

    args = [
        'STAR',
        '--runMode', 'alignReads',
        '--runThreadN', '{:d}'.format(threads),
        '--genomeDir', '{:s}'.format(genome_index),
        '--readFilesIn', '{:s}'.format(sequences_fname2),
    ]
    if not out_dir.endswith('/'):
        out_dir += '/'

    args.extend([
        '--outFileNamePrefix', '{:s}'.format(out_dir),
        '--outSAMprimaryFlag', 'AllBestScore',
        '--outFilterMultimapNmax', '{:d}'.format(multimax),
        '--outFilterMismatchNmax', '{:d}'.format(mismatches),
        '--alignEndsType', 'EndToEnd',    # soft-clipping of starts and ends may produce multiple hits
        '--outSAMtype', 'BAM', 'SortedByCoordinate',
        '--outSAMunmapped', 'Within', 'KeepPairs',
    ])
    if genome_load:
        args.extend([
            '--genomeLoad', 'LoadAndRemove',  # load genome into shared memory
        ])

    if annotation:
        annotation2 = iCount.files.decompress_to_tempfile(annotation, 'starmap')
        args.extend([
            '--sjdbGTFfile', annotation2,
        ])
    else:
        annotation2 = annotation

    try:
        ret_code = 1
        for name, value in _execute(args):
            if name == 'return_code' and isinstance(value, int):
                ret_code = value
                break
            elif name == 'stdout_line':
                LOGGER.info(value.strip())
            elif name == 'stderr_lines':
                for line in value.split('\n'):
                    LOGGER.error(line.strip())
    finally:
        # remove temporary decompressed files
        if reads != sequences_fname2:
            os.remove(sequences_fname2)
        if annotation != annotation2:
            os.remove(annotation2)

    LOGGER.info('Done.')
    return ret_code
