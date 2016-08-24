"""
Map
---

Calls STAR to map sequence reads to reference genome, by passing the
following parameters:

--runThreadN NumberOfThreads
--genomeDir /path/to/genomeDir
--readFilesIn /path/to/read1 [/path/to/read2]
--readFilesCommand gunzip -c
--sjdbGTFfile /path/to/ann.gtf
--outFileNamePrefix /path/to/output/dir/prefix
#--outSAMmapqUnique Integer0to255
--outSAMprimaryFlag AllBestScore
--outFilterMultimapNmax maxHits # default 10
--outFilterMismatchNmax # default 10
"""

import os
import iCount

# description and parameters needed for the analysis
analysis_name = 'map'
analysis_description_short = 'map to genome'
analysis_description = 'Call STAR to map reads to genome index and produce ' \
                       'BAM file.'

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
        'multimax', 'int_range', (50, 1, 300), False,
        'Number of allowed multiple hits.'
    ),
    (
        'mismatches', 'int_range', (2, 0, 20), False,
        'Number of allowed mismatches.'
    ),
]

params_pos = [
    (
        'sequences', 'FASTQ', 'in',
        'Sequencing reads to map to genome.'
    ),
    (
        'genome', 'folder', 'in',
        'Folder with genome index.'
    ),

]


def run(sequences_fname, genomedir, outdir, annotation_fname='', multimax=50,
        mismatches=2, threads=1):
    assert os.path.isdir(outdir)
    assert os.path.isdir(genomedir)
    return iCount.externals.star.map_reads(sequences_fname, genomedir,
                                           outdir,
                                           annotation=annotation_fname,
                                           multimax=multimax,
                                           mismatches=mismatches,
                                           threads=threads)
