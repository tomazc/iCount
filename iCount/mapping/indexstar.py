""".. Line to protect from pydocstyle D205, D400.

Mapping index
-------------

Generate STAR index for mapping based on genome sequence and annotation.

Calls STAR to generate index by passing the following parameters:

--runThreadN NumberOfThreads
--runMode genomeGenerate
--genomeDir /path/to/genomeDir
--genomeFastaFiles /path/to/genome/fasta1 /path/to/genome/fasta2 ...
--sjdbGTFfile /path/to/annotations.gtf
--sjdbOverhang ReadLength-1

"""
# pylint: disable=unused-import
from iCount.externals.star import build_index as run
