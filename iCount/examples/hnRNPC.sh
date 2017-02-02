#!/usr/bin/env bash
set -vx

##############################
# create "common" folders
genome_dir="genomes"
mkdir -p ${genome_dir}

# create analysis output folder
analysis_dir="hnRNPC"
mkdir -p ${analysis_dir}


##############################
# ## general steps
# step 1: download genome sequence
# if parameter --release is omitted, the newest supported release is used
genome="hs84.fa.gz"
if [ ! -f "${genome_dir}/${genome}" ];
then
iCount genome homo_sapiens --release 84 --target_dir ${genome_dir} \
    --target_fname ${genome}
fi
genome="${genome_dir}/${genome}"


# step 2: download annotation
# if parameter --release is omitted, the newest supported release is used
annotation="hs84.gtf.gz"
if [ ! -f "${genome_dir}/${annotation}" ];
then
    iCount annotation homo_sapiens --release 84 --target_dir ${genome_dir} \
        --target_fname ${annotation}
fi
annotation="${genome_dir}/${annotation}"

segmentation="${genome_dir}/hs84.genes.gtf.gz"
if [ ! -f "${segmentation}" ];
then
    iCount genes ${annotation} ${segmentation}
fi


# step 3: index genome
genome_index="star_index/hg19"
if [ ! -d ${genome_index} ];
then
    mkdir -p ${genome_index}
    overhang="100"
    iCount indexstar ${genome} ${genome_index} --annotation ${annotation} \
        --overhang ${overhang} --threads 1
fi


##############################
### analysis-specific steps
pushd ${analysis_dir}

# step 4: download FASTQ file of a hnRNPC iCLIP library
fastq="hnRNPC.fq.gz"
if [ ! -f ${fastq} ];
then
    url="http://icount.fri.uni-lj.si/data/20101116_LUjh03/SLX-2605\
.CRIRUN_501.s_4.sequence.txt.gz"
    wget --no-verbose -O ${fastq} ${url}
fi


# step 5: demultiplex
barcodes="NNNGGTTNN, NNNTTGTNN, NNNCAATNN, NNNACCTNN, NNNGGCGNN"
adapter="AGATCGGAAGAGCGGTTCAG"
mm="1"  # number of mismatches allowed in sample barcode
ext_pref="exp"
ext_dir="extracted"
mkdir -p ${ext_dir}
iCount demultiplex ${fastq} ${adapter} ${barcodes} --mismatches ${mm} --prefix ${ext_pref} \
    --outdir ${ext_dir}


# step 6: map reads to genome reference
multimax="50"
mismatches="2"

for barcode in ${barcodes//,/ }; do
    sequences="${ext_dir}/${ext_pref}_${barcode}.fastq.gz"

    map_dir="mappings/${ext_pref}_${barcode}"
    rm -Rf ${map_dir}
    mkdir -p ${map_dir}

    iCount mapstar ${sequences} "../${genome_index}" ${map_dir} --threads 1 \
        --multimax ${multimax} --mismatches ${mismatches} \
        --annotation "../${annotation}"
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

    iCount xlsites ${bam} ${bed} ${bedm} --mismatches ${mismatches} \
        --group_by ${groupby} --quant ${quant}


    # step 8: perform peaks analysis
    peaks="${peaks_dir}/peaks_${barcode}_${groupby}_${quant}_unique.bed"
    scores="${peaks_dir}/peaks_${barcode}_${groupby}_${quant}_unique.tab"
    iCount peaks "../${segmentation}" ${bed} ${peaks}
        --fout_scores ${scores}
done

popd
