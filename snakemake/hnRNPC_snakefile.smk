#==============================================================================#
#                        iCount Snakemake workflow
#==============================================================================#
# Authors # Igor Ruiz de los Mozos, Charlotte Capitanchik, Tomaz Curk
# Last updated: August 2019


# HOW TO RUN

# Step one: Activate conda environment with Snakemake, iCount and dependencies installed
# Create new environment
# conda env create --name iCount_pipeline --file envs/environment_iCount.yaml
# conda activate iCount_pipeline
# pip install ./iCount/
# Check the install
# iCount

# Step two: To run locally use command:
# snakemake -k -p --snakefile hnRNPC_snakefile.smk --use-conda


# ENVS by Charlotte
# conda create -n iCount_pipeline python=3
# conda install -n iCount_pipeline -c bioconda -c conda-forge snakemake anaconda bedtools bowtie STAR trim-galore sambamba deeptools pip pybedtools samtools=1.8 pytest conda

# ENVS by IGOR on CAMP
# ml Anaconda2/2019.03
# conda create -n iCount_pipeline python=3.6
# conda install -n iCount_pipeline -c bioconda snakemake jinja2 networkx bcftools samtools pysam cutadapt bedtools STAR pip numpy pandas pybedtools numpydoc sphinx matplotlib docutils sphinx_rtd_theme
# conda activate iCount_pipeline
# cd /camp/lab/ulej/working/Igor/Programs/iCount
# pip install ./iCount/
# conda update snakemake

# ENVS2 by IGOR on CAMP
# ml Python/3.6.6-foss-2018b
# ml Anaconda2/2019.03
# conda-env create --file envs/environment_iCount.yaml
# source activate iCount_pipeline2
# cd /camp/lab/ulej/working/Igor/Programs
# pip install ./iCount/


# Remove everyting
# rm -r annotated/ demultiplexed/ genomes/ logs/ mapped/ metrics/ qc/ sig_xlsites/ trimmed/ xlsites/ untitled\ folder*

# To run in a cluster use command:
# sbatch -J iCount_main -o iCount_%A.out -N 1 -t 3-00:00:00 --wrap="snakemake -k -p --snakefile hnRNPC_snakefile.smk --jobs 99 --use-conda --cluster-config cluster_slurmn.yaml --cluster 'sbatch -J {cluster.name} -N {cluster.n} -c {cluster.c} --mem={cluster.memory} -t {cluster.time} -o {cluster.output} -e {cluster.error}'"

# conda activate iCount_pipeline
# mkdir logs
#

# Dry run
# snakemake -k -p -n -r --snakefile hnRNPC_snakefile.smk --use-conda
# Unlock directory
# snakemake --unlock -k -p --snakefile hnRNPC_snakefile.smk

from os.path import join

shell.executable("/bin/bash")

import re
import pandas as pd

#from snakemake.utils import validate



#~~~~~~~~~~~~~~~~~~~~~* Import config file and parameters *~~~~~~~~~~~~~~~~~~~~#
configfile:"config.yaml"
# validate(config, schema="schemas/config.schema.yaml")

samples = pd.read_table(config["samples"]).set_index("5' barcode", drop=False)
#validate(samples, schema="schemas/samples.schema.yaml")

#~~~~~~~~~~~~~~~~~~~~~* Create log folder for cluster run *~~~~~~~~~~~~~~~~~~~~#
logdir = os.path.join(os.getcwd(), config["logdir"])
os.makedirs(logdir, exist_ok=True)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~* Final outputs *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
localrules: all, create_cluster_log

rule all:
    input:
        "demultiplexed/demux_nomatch.fastq.gz",
        "qc/fastqc/raw_fastq_file_fastqc.html",
        "qc/fastqc/raw_fastq_file_fastqc.zip",
        expand("qc/fastqc/{barcode}_fastqc.html", barcode=samples.index),
        expand("qc/fastqc/{barcode}_fastqc.zip", barcode=samples.index),
        expand("trimmed/demux_{barcode}_trimmed.fastq.gz", barcode=samples.index),
        expand("qc/fastqc/{barcode}_trimmed_fastqc.html", barcode=samples.index),
        expand("qc/fastqc/{barcode}_trimmed_fastqc.zip", barcode=samples.index),
        expand("genomes/{genome}/{genome}.fa.gz", genome=samples["mapto"].unique()),
        expand("genomes/{genome}/{genome}.fa.gz.fai", genome=samples["mapto"].unique()),
        expand("genomes/{genome}/{genome}.gtf.gz", genome=samples["mapto"].unique()),
        expand("genomes/{genome}/star_index/", genome=samples["mapto"].unique()),
        expand("mapped/{barcode}/Aligned.sortedByCoord.out.bam", barcode=samples.index),
        expand("xlsites/{barcode}.unique.xl.bed", barcode=samples.index),
        expand("genomes/{genome}/segment/{genome}_segment.gtf", genome=samples["mapto"].unique()),
        expand("sig_xlsites/{barcode}.lowFDR.bed", barcode=samples.index),
        expand("annotated/{barcode}.unique.xl.annotated_sites_biotype.tab", barcode=samples.index),
        expand("annotated/{barcode}.unique.xl.annotated_sites_gene_id.tab", barcode=samples.index),
        expand("annotated/summary_{barcode}/summary_gene.tsv", barcode=samples.index),
        # expand("logs/{barcode}_mapstar.Log.final.out",  barcode=samples.index),


# output:
# genome_fasta = "genomes/{genome}/{genome}.fa.gz",
# genome_index = "genomes/{genome}/{genome}.fa.gz.fai",
# gtf = "genomes/{genome}/{genome}.gtf.gz"



#ruleorder: demultiplex > quality_trim > fastqc_trimmed > download_annotation > index_genome > map_reads


#==============================================================================#
#                       Demultiplex
#==============================================================================#




#### Include --mismatches to barcode demultiplex
#### Check if one of the barcodes OR nomatch is not created/found
#### Include 3' barcode demultiplex
rule demultiplex:
    input:
        fastq_file=config['raw_fastq_file']
    output:
        expand("demultiplexed/demux_{barcode}.fastq.gz", barcode=samples.index),
        "demultiplexed/demux_nomatch.fastq.gz"
    params:
        adapter3=config['adapter3'],
        # cluster='-N 1 -c 1 -J demultiplex --mem=32G -t 16:00:00',
        all_5barcodes = samples["5' barcode"].unique().tolist()
    log:
        metrics="metrics/demultiplex_metrics.txt"
    shell:
        """
        iCount demultiplex {input.fastq_file} {params.adapter3} {params.all_5barcodes} --out_dir demultiplexed -M {log.metrics}
        """


#==============================================================================#
#                       Read quality trimming and QC
#==============================================================================#
# def construct name of raw fastq config['raw_fastq_file']

rule fastqc_raw:
    input:
        config['raw_fastq_file']
    output:
        html="qc/fastqc/raw_fastq_file_fastqc.html",
        zip="qc/fastqc/raw_fastq_file_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    # params:
    #     cluster='-N 1 -c 5 -J raw_fastqc --mem=30G -t 4:00:00'
    log:
        "logs/fastqc/raw_fastq_file_fastqc.log"
    wrapper:
        "0.36.0/bio/fastqc"


rule fastqc:
    input:
        "demultiplexed/demux_{barcode}.fastq.gz"
    output:
        html="qc/fastqc/{barcode}_fastqc.html",
        zip="qc/fastqc/{barcode}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    # params:
    #     cluster= '-N 1 -c 5 -J fastqc --mem=30G -t 4:00:00'
    log:
        "logs/fastqc/{barcode}_fastqc.txt"
    wrapper:
        "0.36.0/bio/fastqc"


#Ask Tomaz why metrics is equal to 0 and the log only contains the cutadapt command
rule quality_trim:
    input:
        "demultiplexed/demux_{barcode}.fastq.gz"
    output:
        trimmed_reads="trimmed/demux_{barcode}_trimmed.fastq.gz",
        metrics="metrics/{barcode}_trimmed.txt"
    params:
        qual_trim=config['qual_trim'],
        minimum_length=config['minimum_length'],
        adapter=config['adapter3'],
        # cluster='-N 1 -c 5 -J quality_trim --mem=30G -t 16:00:00'
    log:
        "logs/trimmed/{barcode}_trimmed.txt"
    shell:
        """
        iCount cutadapt --qual_trim {params.qual_trim} --minimum_length {params.minimum_length} --file_log 2 --file_logpath {log} --results_file {output.metrics} {input} {output.trimmed_reads} {params.adapter}
        """


rule fastqc_trimmed:
    input:
        "trimmed/demux_{barcode}_trimmed.fastq.gz"
    output:
        html="qc/fastqc/{barcode}_trimmed_fastqc.html",
        zip="qc/fastqc/{barcode}_trimmed_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    # params:
    #     cluster='-N 1 -c 5 -J fastqc_trimmed --mem=30G -t 4:00:00'
    log:
        "logs/fastqc/{barcode}_trimmed_fastqc.log"
    wrapper:
        "0.36.0/bio/fastqc"


#==============================================================================#
#                       Download annotation and index genome
#==============================================================================#


# if the genome is not there it should complain and fail with hint to download annotation. This to validation of tabular sample file!!
# Only once there was an error running on the cluster (broken connection? NCBI issue?) "RuleException: CalledProcessError" it could be related to remote files download "https://snakemake.readthedocs.io/en/stable/snakefiles/remote_files.html"

rule download_annotation:
    output:
        genome_fasta="genomes/{genome}/{genome}.fa.gz",
        genome_index="genomes/{genome}/{genome}.fa.gz.fai",
        gtf="genomes/{genome}/{genome}.gtf.gz"
    params:
        # cluster='-N 1 -c 1 -J download_annotation --mem=8G -t 4:00:00',
        release=config['release']
    shell:
        """
        iCount genome --genome {output.genome_fasta} --source ensembl {wildcards.genome} {params.release} --chromosomes MT 19
        iCount annotation --annotation {output.gtf} --source ensembl {wildcards.genome} {params.release}
        """



rule index_genome:
    input:
        genome_fasta="genomes/{genome}/{genome}.fa.gz",
        gtf="genomes/{genome}/{genome}.gtf.gz"
    threads:
        8
    params:
        overhang=config['ref']['overhang'],
        # cluster='-J star_index -t 16:00:00 -n 1 -c 1'
        #--mem-per-cpu 16GB -c 8
    output:
        directory("genomes/{genome}/star_index/")
    shell:
        """
        iCount indexstar --overhang {params.overhang} --annotation {input.gtf} \
        --threads {threads} --genome_sasparsed 2 {input.genome_fasta} {output}
        """

#==============================================================================#
#                       Map reads
#==============================================================================#


def get_gtf_path(wildcards):
    return ("genomes/{0}/{0}.gtf.gz".format(samples.loc[wildcards.barcode, "mapto"]))

def get_star_index_path(wildcards):
    return ("genomes/{}/star_index/".format(samples.loc[wildcards.barcode, "mapto"]))

def get_segment_path(wildcards):
    return ("genomes/{0}/segment/{0}_segment.gtf".format(samples.loc[wildcards.barcode, "mapto"]))

def get_templates_dir(wildcards):
    return ("genomes/{0}/segment/".format(samples.loc[wildcards.barcode, "mapto"]))



# Include bam.bai
rule map_reads:
    input:
        trimmed_reads="trimmed/demux_{barcode}_trimmed.fastq.gz",
        star_index = directory(get_star_index_path),
        gtf = get_gtf_path,
    output:
        "mapped/{barcode}/Aligned.sortedByCoord.out.bam"
    params:
        outdir=directory("mapped/{barcode}/"),
        multimax=config['multimax'],
        # cluster='-J mapstar  -t 16:00:00 -n 1 --mem-per-cpu 8GB -c 8'
    log:
        "logs/mapstar/{barcode}_mapstar.log"
    shell:
        """
        iCount mapstar --annotation {input.gtf} --multimax {params.multimax} \
        {input.trimmed_reads} {input.star_index} {params.outdir}
        """


#==============================================================================#
#                       Create cross link sites and peaks
#==============================================================================#


rule xlsites:
    input:
        "mapped/{barcode}/Aligned.sortedByCoord.out.bam"
    output:
        unique_bed="xlsites/{barcode}.unique.xl.bed",
        multimapped_bed="xlsites/{barcode}.multimapped.xl.bed",
        skipped_bam="xlsites/{barcode}.skipped.xl.bam",
    # params:
    #     cluster='-N 1 -c 1 --mem=30G -t 16:00:00'
    benchmark:
        "benchmarks/{barcode}.xlsites.benchmark.tab"
        # repeat("benchmarks/{barcode}.xlsites.benchmark.tab", 3)
    log:
        "logs/xlsites/{barcode}.xlsites.log"
    shell:
        """
        iCount xlsites {input} {output.unique_bed} {output.multimapped_bed} {output.skipped_bam} -M {log}
        """

# rule bedgraph: ??


rule segment:
    input:
        gtf="genomes/{genome}/{genome}.gtf.gz",
        genome_fai="genomes/{genome}/{genome}.fa.gz.fai"
    output:
        "genomes/{genome}/segment/{genome}_segment.gtf"
    # params:
    #     cluster='-N 1 -c 1 --mem=30G -t 16:00:00'
    shell:
        """
        iCount segment {input.gtf} {output} {input.genome_fai}
        """


rule sig_xlsites:
    input:
        xlsites="xlsites/{barcode}.unique.xl.bed",
        segment_file=get_segment_path
    output:
        sigxls="sig_xlsites/{barcode}.lowFDR.bed",
        scores="sig_xlsites/{barcode}.scores.tsv"
    # params:
    #     cluster='-N 1 -c 1 --mem=30G -t 16:00:00'
    benchmark:
        "benchmarks/{barcode}.sig_xlsites.benchmark.txt"
        # repeat("benchmarks/{barcode}.sig_xlsites.benchmark.tab", 3)
    shell:
        """
        iCount peaks {input.segment_file} {input.xlsites} {output.sigxls} --scores {output.scores}
        """

# rule bedgraph: ??


#==============================================================================#
#                       Annotate cross link sites and summaries
#==============================================================================#


rule annotate_xlsites:
    input:
        xlsites = "xlsites/{barcode}.unique.xl.bed"
    output:
        biotype="annotated/{barcode}.unique.xl.annotated_sites_biotype.tab",
        gene_id="annotated/{barcode}.unique.xl.annotated_sites_gene_id.tab"
    params:
        templates_dir = get_templates_dir,
        segment = get_segment_path,
        out_dir = "annotated/",
        # cluster='-N 1 -c 1 --mem=30G -t 16:00:00'
    shell:
        """
        iCount annotate --subtype biotype {params.segment} {input.xlsites} {output.biotype}
        iCount annotate --subtype gene_id {params.segment} {input.xlsites} {output.gene_id}
        """


# Summary on xlsites or peaks??
rule summary:
    input:
        xlsites="xlsites/{barcode}.unique.xl.bed"
    output:
        "annotated/summary_{barcode}/summary_gene.tsv"
    # conda:
    #     "envs/summary.yaml"
    params:
        templates_dir=get_templates_dir,
        gtf = get_gtf_path,
        segment = get_segment_path,
        out_dir="annotated/summary_{barcode}/",
        # cluster='-N 1 -c 1 --mem=30G -t 16:00:00'
    shell:
        """
        iCount summary --templates_dir {params.templates_dir} {params.segment} {input.xlsites} {params.out_dir}
        """

#==============================================================================#
#                       Group analysis
#==============================================================================#


# rule group:


#==============================================================================#
#                       RNAmaps & Kmers
#==============================================================================#




# rule RNAmaps:
# rule kmers:
# rule clean:

#==============================================================================#
#                       Metrics and experiment QC
#==============================================================================#

# multiQC implementation
# total reads
# demultiplex %
# total xlistes % of PCR duplication
# total peaks
# biotypes per sample
# Lits of more abundant binders
# More...?