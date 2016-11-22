""".. Line to protect from pydocstyle D205, D400.

Mapping reads with STAR
-----------------------

Map reads to genome index and produce BAM file by using STAR.

Call STAR by passing the following parameters:

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
# pylint: disable=unused-import
from iCount.externals.star import map_reads as run
