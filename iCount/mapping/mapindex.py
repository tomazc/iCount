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
from iCount.externals.star import build_index as run
