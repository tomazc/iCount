"""
Call STAR to map reads to genome index and produce BAM file.

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

from iCount.externals.star import map_reads as run
