#!/usr/bin/env bash
set -vx

##############################
# create "common" folders
genome_dir="genomes"
mkdir -p ${genome_dir}

# create analysis output folder
analysis_dir="hnRNPCred"
mkdir -p ${analysis_dir}


##############################
# ## general steps
# step 1: download genome sequence
genome="hs84red.fa.gz"
if [ ! -f "${genome_dir}/${genome}" ];
then
iCount sequence -release 84 -species homo_sapiens -chromosomes "21,MT" \
                -target_dir ${genome_dir} -target_fname ${genome}

mv "${genome_dir}/hs84red.chr21_MT.fa.gz" "${genome_dir}/${genome}"
mv "${genome_dir}/hs84red.chr21_MT.fa.gz.chrom_length.txt" \
   "${genome_dir}/hs84red.fa.gz.chrom_length.txt"
fi
genome="${genome_dir}/${genome}"


# step 2: download annotation
annotation="hs84.gtf.gz"
if [ ! -f "${genome_dir}/${annotation}" ];
then
    iCount annotation -release 84 -species homo_sapiens \
                      -target_dir ${genome_dir} -target_fname ${annotation}
fi
annotation="${genome_dir}/${annotation}"

genes_annotation="${genome_dir}/hs84.genes.bed.gz"
if [ ! -f "${genes_annotation}" ];
then
    iCount segment ${annotation} ${genes_annotation}
fi


# step 3: index genome
genome_index="star_index/hg19red"
if [ ! -d ${genome_index} ];
then
    mkdir -p ${genome_index}
    overhang="100"
    iCount mapindex -threads 1 -outdir ${genome_index} \
                    -annotation ${annotation} -overhang ${overhang} ${genome}
fi


##############################
### analysis-specific steps
pushd ${analysis_dir}

# step 4: download FASTQ file of a hnRNPC iCLIP library
fastq="hnRNPC.fq.gz"
if [ ! -f ${fastq} ];
then
    url="http://icount.fri.uni-lj.si/data/20101116_LUjh03/SLX-2605\
.CRIRUN_501.s_4.sequence.reduced.txt.gz"
    wget --no-verbose -O ${fastq} ${url}
fi


# step 5: demultiplex
barcodes="NNNGGTTNN,NNNTTGTNN,NNNCAATNN,NNNACCTNN,NNNGGCGNN"
adapter="AGATCGGAAGAGCGGTTCAG"
mm="1"  # number of mismatches allowed in sample barcode
ext_pref="exp"
ext_dir="extracted"
mkdir -p ${ext_dir}
iCount demultiplex -barcodes ${barcodes} -adapter ${adapter} \
                   -mismatches ${mm} -prefix ${ext_pref} -outdir ${ext_dir} \
                   ${fastq}


# step 6: map reads to genome reference
multimax="50"
mismatches="2"

for barcode in ${barcodes//,/ }; do
    sequences="${ext_dir}/${ext_pref}_${barcode}.fastq.gz"

    map_dir="mappings/${ext_pref}_${barcode}"
    rm -Rf ${map_dir}
    mkdir -p ${map_dir}

    iCount map -threads 1 -outdir ${map_dir} -multimax ${multimax} \
               -mismatches ${mismatches} ${sequences} "../${genome_index}"
done


# step 7: identify and quantify cross-linked sites
groupby="start"
quant="cDNA"
mismatches="2"

bed_dir="beds"
rm -Rf ${bed_dir}
mkdir -p ${bed_dir}

peaks_dir="peaks"
rm -Rf ${peaks_dir}
mkdir -p ${peaks_dir}

for barcode in ${barcodes//,/ }; do
    map_dir="mappings/${ext_pref}_${barcode}"
    bam="${map_dir}/Aligned.sortedByCoord.out.bam"

    bed="${bed_dir}/hits_${barcode}_${groupby}_${quant}_unique.bed"
    bedm="${bed_dir}/hits_${barcode}_${groupby}_${quant}_multi.bed"

    iCount xlsites ${bam} ${bed} ${bedm} -mismatches ${mismatches} \
                   -groupby ${groupby} -quant ${quant}


    # step 8: perform peaks analysis
    peaks="${peaks_dir}/peaks_${barcode}_${groupby}_${quant}_unique.bed"
    scores="${peaks_dir}/peaks_${barcode}_${groupby}_${quant}_unique.tab"
    iCount peaks "../${genes_annotation}" ${bed} ${peaks} ${scores} \
                 -verbose no
done

popd
