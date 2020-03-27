#!/usr/bin/env bash
# Create genome, annotation and read data for an iCLIP test data set
# Charlotte Capitanchik 13.03.2020

chr="chr6"
start="34,000,000"
end="35,000,000"
start_annot="34,100,000"
end_annot="34,900,000"
bam_dir="/camp/lab/luscomben/working/charlotte/SM/hg38_Gencode29_BPs_LaBranchoR/PRP8_rupert/iMaps_bams/"
genome_fq="/camp/lab/luscomben/working/charlotte/spliceosome_iCLIP/permanent_files/GRCh38.primary_assembly.genome.fa"
genome_annot="/camp/lab/luscomben/working/charlotte/hg38_regions/gencode.v30.primary_assembly.annotation.gtf"

# ~~ To get the reads ~~ #
bam_files="$(ls ${bam_dir} | grep bam$)"

for bam in $bam_files; do
    filename="$(basename $bam | sed s'/\.bam//g')"
    newbam=${filename}_${chr}_${start//,}_${end//,}.bam
    newfastq=${filename}_${chr}_${start//,}_${end//,}.fq
    samtools view -h ${bam_dir}${bam} ${chr}:${start}-${end} > $newbam
    echo "made $newbam"
    bedtools bamtofastq -i $newbam -fq $newfastq
    echo "made $newfastq"
    gzip $newfastq
    echo "zipped $newfastq"
done

# ~~ To get the genome fastq ~~ #
newgenome="$(basename ${genome_fq} | sed s'/\.fa//g')"
samtools faidx ${genome_fq} ${chr}:${start}-${end} > ${newgenome}_${chr}_${start//,}_${end//,}.fa

# ~~ To get the genome annotation ~~ #
# First make bed file of region
echo -e "${chr}\t${start_annot//,}\t${end_annot//,}" > region.bed
echo -e "${chr}\t${start_annot//,}\t${end_annot//,}" >> region.bed
newannot="$(basename ${genome_annot} | sed s'/\.gtf//g')"
bedtools intersect -a ${genome_annot} -b region.bed -u > ${newannot}_${chr}_${start//,}_${end//,}.gtf