#!/usr/bin/env bash
set -vx

mkdir tutorial_example
cd tutorial_example

iCount releases

iCount species -r 84

iCount genome homo_sapiens -r 84 --chromosomes 21 MT

iCount annotation homo_sapiens -r 84

mkdir hs84
iCount indexstar homo_sapiens.84.chr21_MT.fa.gz hs84 --annotation homo_sapiens.84.gtf.gz

wget http://icount.fri.uni-lj.si/data/20101116_LUjh03/\
SLX-2605.CRIRUN_501.s_4.sequence.txt.gz -O hnRNPC.fq.gz

mkdir demultiplexed
iCount demultiplex hnRNPC.fq.gz AGATCGGAAGAGCGGTTCAG NNNGGTTNN NNNTTGTNN \
NNNCAATNN NNNACCTNN NNNGGCGNN --out_dir "demultiplexed"

ls -lh demultiplexed

mkdir mapping_NNNACCTNN
iCount map demultiplexed/demux_NNNACCTNN.fastq.gz hs84 mapping_NNNACCTNN \
--annotation homo_sapiens.84.gtf.gz

ls -lh mapping_NNNACCTNN

iCount xlsites mapping_NNNACCTNN/Aligned.sortedByCoord.out.bam \
NNNACCTNN_cDNA_unique.bed  NNNACCTNN_cDNA_multiple.bed NNNACCTNN_cDNA_skipped.bam \
--group_by start --quant cDNA

iCount xlsites mapping_NNNACCTNN/Aligned.sortedByCoord.out.bam \
NNNACCTNN_reads_unique.bed  NNNACCTNN_reads_multiple.bed NNNACCTNN_reads_skipped.bam \
--group_by start --quant reads

iCount segment homo_sapiens.84.gtf.gz hs84seg.gtf.gz \
homo_sapiens.84.chr21_MT.fa.gz.fai

iCount peaks hs84seg.gtf.gz NNNACCTNN_cDNA_unique.bed peaks.bed \
--scores scores.tsv
