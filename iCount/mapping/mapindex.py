"""
Call STAR to generate genome index, which is used for mapping.

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


def run(genome_fname, outdir, annotation_fname='', overhang=100, threads=1):
    """
    Call STAR to generate genome index, which is used for mapping.

    Parameters
    ----------
    genome_fname : str
        Genome sequence to index.
    outdir : str
        Output folder, where to store genome index.
    annotation_fname : str
        Annotation that defines splice junctions.
    overhang : int
        Sequence length around annotated junctions to be used by STAR when
        constructing splice junction database.
    threads : int
        Number of threads that STAR can use for generating index.

    Returns
    -------
    int
        Star return code.

    """
    assert os.path.isdir(outdir)
    return iCount.externals.star.build_index(
        genome_fname,
        outdir,
        annotation=annotation_fname,
        overhang=overhang,
        threads=threads,
    )
