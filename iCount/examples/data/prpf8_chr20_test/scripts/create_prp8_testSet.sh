# Create genome, annotation and read data for an iCLIP test data set
# Charlotte Capitanchik 13.03.2020
chr="chr20"
bam_dir="/camp/lab/luscomben/home/users/capitac/SM/hg38_Gencode29_BPs_LaBranchoR/PRP8_rupert/iMaps_bams/"
genome_fq="/camp/lab/luscomben/home/users/capitac/spliceosome_iCLIP/permanent_files/GRCh38.primary_assembly.genome.fa"
genome_annot="/camp/lab/luscomben/home/users/capitac/hg38_regions/gencode.v30.primary_assembly.annotation.gtf"
# ~~ To get the reads ~~ #
bam_files="$(ls ${bam_dir} | grep bam$)"
for bam in $bam_files; do
    filename="$(basename $bam | sed s'/\.bam//g')"
    newbam=${filename}_${chr}.bam
    newfastq=${filename}_${chr}.fq
    samtools view -h ${bam_dir}${bam} ${chr} > $newbam
    echo "made $newbam"
    bedtools bamtofastq -i $newbam -fq $newfastq
    echo "made $newfastq"
    gzip $newfastq
    echo "zipped $newfastq"
done
# ~~ To get the genome fastq ~~ #
newgenome="$(basename ${genome_fq} | sed s'/\.fa//g')"
samtools faidx ${genome_fq} ${chr} > ${newgenome}_${chr}.fa