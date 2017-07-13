#!/usr/bin/env bash
set -vx

mkdir tutorial_example
cd tutorial_example

iCount releases

iCount species -r 88

iCount genome homo_sapiens -r 88 --chromosomes 21 MT

iCount annotation homo_sapiens -r 88

mkdir hs88
iCount indexstar homo_sapiens.88.chr21_MT.fa.gz hs88 --annotation homo_sapiens.88.gtf.gz

# the whole data set [880 MB] is available here:
#wget http://icount.fri.uni-lj.si/data/20101116_LUjh03/\
#SLX-2605.CRIRUN_501.s_4.sequence.txt.gz -O hnRNPC.fq.gz

# in this tutorial, we are using a subset of reads [23 MB]
wget http://icount.fri.uni-lj.si/data/20101116_LUjh03/SLX-2605\
.CRIRUN_501.s_4.sequence.reduced.txt.gz -O hnRNPC.fq.gz

mkdir demultiplexed
iCount demultiplex hnRNPC.fq.gz AGATCGGAAGAGCGGTTCAG NNNGGTTNN NNNTTGTNN \
NNNCAATNN NNNACCTNN NNNGGCGNN --out_dir demultiplexed

ls -lh demultiplexed

mkdir mapping_NNNGGCGNN
iCount mapstar demultiplexed/demux_NNNGGCGNN.fastq.gz hs88 mapping_NNNGGCGNN \
--annotation homo_sapiens.88.gtf.gz

ls -lh mapping_NNNGGCGNN

iCount xlsites mapping_NNNGGCGNN/Aligned.sortedByCoord.out.bam \
NNNGGCGNN_cDNA_unique.bed  NNNGGCGNN_cDNA_multiple.bed NNNGGCGNN_cDNA_skipped.bam \
--group_by start --quant cDNA

iCount xlsites mapping_NNNGGCGNN/Aligned.sortedByCoord.out.bam \
NNNGGCGNN_reads_unique.bed  NNNGGCGNN_reads_multiple.bed NNNGGCGNN_reads_skipped.bam \
--group_by start --quant reads

iCount segment homo_sapiens.88.gtf.gz hs88seg.gtf.gz \
homo_sapiens.88.chr21_MT.fa.gz.fai

iCount peaks hs88seg.gtf.gz NNNGGCGNN_cDNA_unique.bed peaks.bed \
--scores scores.tsv

iCount clusters peaks.bed clusters.bed

iCount annotate hs88seg.gtf.gz NNNGGCGNN_cDNA_unique.bed annotated_sites_biotype.tab
iCount annotate --subtype gene_id hs88seg.gtf.gz NNNGGCGNN_cDNA_unique.bed \
annotated_sites_genes.tab

iCount summary hs88seg.gtf.gz NNNGGCGNN_cDNA_unique.bed summary.tab \
homo_sapiens.88.chr21_MT.fa.gz.fai
