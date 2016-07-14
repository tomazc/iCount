[![Build Status](https://travis-ci.com/tomazc/iCount.svg?token=MxKtDvsXZMsCDvfFpmd6&branch=master)](https://travis-ci.com/tomazc/iCount)
# iCount processing of iCLIP protein-RNA interaction data

This is the command-line interface (CLI) to the iCount pipeline.
The pipeline provides all the commands needed to process FASTQ iCLIP
data and generate:
 
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
