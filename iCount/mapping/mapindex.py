"""
Map index
---------

Calls STAR to generate index based on genome sequence and annotation,
by passing the following parameters:

--runThreadN NumberOfThreads
--runMode genomeGenerate
--genomeDir /path/to/genomeDir
--genomeFastaFiles /path/to/genome/fasta1 /path/to/genome/fasta2 ...
--sjdbGTFfile /path/to/annotations.gtf
--sjdbOverhang ReadLength-1

"""

import os
import iCount

# description and parameters needed for the analysis
analysis_name = 'mapindex'
analysis_description_short = 'generate genome index'
analysis_description = 'Call STAR to generate genome index, which is used ' \
                       'for mapping.'

params_opt = [
    (
        'threads', 'int_range', (1, 1, 32), False,
        'Number of threads that STAR can use for generating index.'
    ),
    (
        'outdir', 'string', 'gdir', False,
        'Output folder, where to store genome index.'
    ),
    (
        'annotation', 'string', '', False,
        'Annotation that defines splice junctions.'
    ),
    (
        'overhang', 'int_range', (100, 30, 300), False,
        'Sequence length around annotated junctions to be used by STAR when '
        'constructing splice junction database.'
    ),
]

params_pos = [
    (
        'genome', 'FASTA', 'in', 1,
        'Genome sequence to index.'
    ),
]


def run(genome_fname, outdir, annotation_fname='', overhang=100, threads=1):
    assert os.path.isdir(outdir)
    return iCount.externals.star.build_index(genome_fname, outdir,
                                             annotation=annotation_fname,
                                             overhang=overhang,
                                             threads=threads)
