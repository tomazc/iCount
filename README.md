[![Build Status](https://travis-ci.com/tomazc/iCount.svg?token=MxKtDvsXZMsCDvfFpmd6&branch=master)](https://travis-ci.com/tomazc/iCount)
[![codecov](https://codecov.io/gh/tomazc/iCount/branch/master/graph/badge.svg?token=JhUJ66rnJ3)](https://codecov.io/gh/tomazc/iCount)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/bb21b3cc5fcd420c885ed12bf8393065)](https://www.codacy.com?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=tomazc/iCount&amp;utm_campaign=Badge_Grade)

# iCount processing of iCLIP protein-RNA interaction data

iCount is a Python module and associated command-line interface (CLI),
which provides all the commands needed to process iCLIP data and 
generate:
 
+ Demultiplexed and adapter-trimmed FASTQ files.
+ BAM files with mapped iCLIP reads.
+ bedGraphs with information on identified protein-RNA cross-linked sites.
+ bedGraphs with statistically significant cross-linked sites (peak finding).
+ bedGraphs with clusters of significant cross-linkeds sites (cluster finding).
+ Grouping of individual replicate experiments.
+ RNAmap generation showing the positional distribution of cross-linked sites relative to genomic landmarks.
+ kmer enrichment analysis.
+ ...

Let's start with a simple example...
