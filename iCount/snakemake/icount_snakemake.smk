#==============================================================================#
#                        iCount Snakemake workflow
#==============================================================================#
# Authors # Igor Ruiz de los Mozos, Charlotte Capitanchik, Tomaz Curk
# Last updated: September 2020


# Install Locally
#================
# Move to your software folder
# git clone https://github.com/tomazc/iCount.git
# vi iCount/iCount/__init__.py
# Change the following line:
# TMP_ROOT = os.environ.get('ICOUNT_TMP_ROOT', 'a directory you have access to')
# TMP_ROOT = os.environ.get('ICOUNT_TMP_ROOT', '/Users/mozosi/Programs/temp_iCount')

# Check the install

# Step one: Activate conda environment with Snakemake, iCount and dependencies installed
# Create new environment
# conda env create --name iCount_pipeline2 --file envs/environment_iCount.yaml
# conda activate iCount_pipeline2
# pip install ./iCount/
# Check the install
# iCount

# Run Locally
#================

# Step two: To run locally use command:
# snakemake -k -p --snakefile demultiplex_snakefile.smk --use-conda
# snakemake -k -p --cores 4 --snakefile demultiplex_snakefile.smk --use-conda
# snakemake -k -p --cores 4 --snakefile '/Users/mozosi/Dropbox (UCL-MN Team)/GitHub/iCount/iCount/snakemake/icount_snakemake.smk' --use-conda --configfile config_synthetic.yaml
# dag workflow
# snakemake --snakefile demultiplex_snakefile.smk --use-conda --dag 2> /dev/null | dot -T png > workflow_bysample.png
# snakemake --snakefile demultiplex_snakefile.smk --use-conda --rulegraph 2> /dev/null | dot -T png > workflow.png

# Install Cluster
#================

# ENVS on CAMP
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

# Run Cluster
#================
# To run in a cluster use command:
# sbatch -J iCount_main -o iCount_%A.out -N 1 -t 3-00:00:00 --wrap="snakemake -k -p --snakefile demultiplex_snakefile.smk --jobs 99 --use-conda --cluster-config envs/cluster_slurmn.yaml --cluster 'sbatch -J {cluster.name} -N {cluster.n} -c {cluster.c} --mem={cluster.memory} -t {cluster.time} -o {cluster.output} -e {cluster.error}'"

# conda activate iCount_pipeline
# mkdir logs
#

# Dry run
# snakemake -k -p -n -r --snakefile iCount_snakefile.smk --use-conda
# Unlock directory
# snakemake --unlock -k -p --snakefile iCount_snakefile.smk


# LOGGER = logging.getLogger(__name__)
from os.path import join

shell.executable("/bin/bash")

import re
import pandas as pd
import numpy as np
import gzip
import os
import shutil
import tempfile
import pysam
import yaml

# Validate config file!!!
from snakemake.utils import validate

#~~~~~~~~~~~~~~~~~~~~~* Import config file and samples annotation *~~~~~~~~~~~~~~~~~~~~#
# Config file can be specified here or on the snakemake call (--configfile config_synthetic.yaml)
# configfile:"config_synthetic.yaml"
# Note that the path to the configfile can be absolute or relative to the directory
# where the Snakefile is run.

# Validate config file
validate(config, schema="schemas/config.schema.yaml")

# Import sample file and validate samples
samples = pd.read_table(config["samples"]).set_index("barcode_5", drop=False)
#samples = pd.read_table(config["samples"])
validate(samples, schema="schemas/samples.schema.yaml")

# Validate adapter and samples integrity
if len(samples["adapter_3"].unique().tolist()) > 1:
    sys.exit("iCount pipeline only accepts a unique 3' adapter")

if len(samples.index) != len(samples["sample_name"].unique().tolist()):
    sys.exit("iCount pipeline only accepts a unique sample names")

# Merge 5'barcode and 3'barcode to create a table index (full barcode)
cols = ['barcode_5', 'barcode_3']
samples["full_barcode"] = samples[cols].apply(lambda x: '_'.join(x.dropna()), axis=1)
samples=samples.set_index(["full_barcode"], drop = False)


# Print project
print("Procesing project:", config['project'], "\n")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~* Final outputs *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

##### load rules #####
include: "rules/common.smk"
include: "rules/demultiplex.smk"
include: "rules/qc.smk"
include: "rules/prepare_genome.smk"
include: "rules/map_reads.smk"
include: "rules/bedgraph_UCSC.smk"
include: "rules/xlsites.smk"
include: "rules/sig_xlsites.smk"
include: "rules/clusters.smk"
include: "rules/group.smk"


##### target rules #####

localrules: all

rule all:
    input:
        "demultiplexed/demux_nomatch5.fastq.gz",

        all_input,
        all_xlsites,
        all_group

# old rule all not used - remove
# "demultiplexed/demux_nomatch5.fastq.gz",
# expand("{genomes_path}/{genome}/{genome}.fa.gz", genome=samples["mapto"].unique(), genomes_path=config['genomes_path']),
# expand("{genomes_path}/{genome}/{genome}.fa.gz.fai", genome=samples["mapto"].unique(), genomes_path=config['genomes_path']),
# expand("{genomes_path}/{genome}/{genome}.gtf.gz", genome=samples["mapto"].unique(), genomes_path=config['genomes_path']),
# expand("{genomes_path}/{genome}/star_index/", genome=samples["mapto"].unique(), genomes_path=config['genomes_path']),
# expand("{genomes_path}/{genome}/segment/{genome}_segment.gtf", genome=samples["mapto"].unique(), genomes_path=config['genomes_path']),
# expand("{genomes_path}/{genome}/segment/landmarks.bed.gz", genome=samples["mapto"].unique(), genomes_path=config['genomes_path']),
#
# expand("qc/fastqc/raw_fastq_file_fastqc.html"),
# expand("qc/fastqc/raw_fastq_file_fastqc.zip"),
# expand("{project}/qc/fastqc/{barcode}_fastqc.html", project=config['project'], barcode=samples.index),
# expand("{project}/qc/fastqc/{barcode}_fastqc.zip", project=config['project'], barcode=samples.index, ),
#
# expand("{project}/trimmed/demux_{barcode}_trimmed.fastq.gz", project=config['project'], barcode=samples.index),
# expand("{project}/qc/fastqc/{barcode}_trimmed_fastqc.html", project=config['project'], barcode=samples.index),
# expand("{project}/qc/fastqc/{barcode}_trimmed_fastqc.zip", project=config['project'], barcode=samples.index),
# expand("{project}/mapped/{barcode}/Aligned.sortedByCoord.out.bam", project=config['project'], barcode=samples.index),

# expand("{project}/xlsites/{barcode}/{barcode}.unique.xl.bed", project=config['project'], barcode=samples.index),
# expand("{project}/xlsites/{barcode}/{barcode}.unique.xl.summary_gene.tsv", project=config['project'],
#        barcode=samples.index),
# expand("{project}/xlsites/{barcode}/{barcode}.unique.xl.annotated_sites_biotype.tab", project=config['project'],
#        barcode=samples.index),
# expand("{project}/xlsites/{barcode}/{barcode}.unique.xl.annotated_sites_gene_id.tab", project=config['project'],
#        barcode=samples.index),
# expand("{project}/xlsites/{barcode}/{barcode}.unique.xl.bedgraph", project=config['project'], barcode=samples.index),
# # expand("{project}/xlsites/{barcode}/{barcode}.unique.xl.UCSC.bedgraph", project=config['project'], barcode=samples.index),
# expand("{project}/xlsites/{barcode}/rnamaps/", project=config['project'], barcode=samples.index),
#
# expand("{project}/sig_xlsites/{barcode}/{barcode}.sig_sites.bed", project=config['project'], barcode=samples.index),
# expand("{project}/sig_xlsites/{barcode}/{barcode}.sig_sites.summary_gene.tsv", project=config['project'],
#        barcode=samples.index),
# expand("{project}/sig_xlsites/{barcode}/{barcode}.sig_sites.annotated.biotype.tab", project=config['project'],
#        barcode=samples.index),
# expand("{project}/sig_xlsites/{barcode}/{barcode}.sig_sites.annotated.gene_id.tab", project=config['project'],
#        barcode=samples.index),
#
# expand("{project}/clusters/{barcode}/{barcode}.clusters.bed", project=config['project'], barcode=samples.index),
# expand("{project}/groups/{group}/xlsites/{group}.group.unique.xl.bed", project=config['project'],
#        group=samples["group"].dropna().unique()),
# expand("{project}/groups/{group}/xlsites/{group}.group.unique.xl.annotated_group_biotype.tab",
#        project=config['project'], group=samples["group"].dropna().unique()),
# expand("{project}/groups/{group}/xlsites/{group}.group.unique.xl.annotated_group_gene_id.tab",
#        project=config['project'], group=samples["group"].dropna().unique()),
# expand("{project}/groups/{group}/xlsites/{group}.group.unique.xl.summary_gene.tsv", project=config['project'],
#        group=samples["group"].dropna().unique()),
# # expand("{project}/groups/{group}/xlsites/{group}.group.unique.xl.summary_type.tsv", project=config['project'], group=samples["group"].dropna().unique()),
# # expand("{project}/groups/{group}/xlsites/{group}.group.unique.xl.summary_subtype.tsv", project=config['project'], group=samples["group"].dropna().unique()),
# expand("{project}/groups/{group}/xlsites/{group}.group.unique.xl.bedgraph", project=config['project'],
#        group=samples["group"].dropna().unique()),
# expand("{project}/groups/{group}/rnamaps/", project=config['project'], group=samples["group"].dropna().unique()),
#
# # expand("{project}/groups/{group}/sig_xlsites/{group}.group.sig_sites.bed", project=config['project'], group=samples["group"].dropna().unique()),
# # expand("{project}/groups/{group}/sig_xlsites/{group}.group.sig_sites.annotated.biotype.tab", project=config['project'], group=samples["group"].dropna().unique()),
# # expand("{project}/groups/{group}/sig_xlsites/{group}.group.sig_sites.annotated.gene_id.tab", project=config['project'], group=samples["group"].dropna().unique()),
# # expand("{project}/groups/{group}/sig_xlsites/{group}.group.sig_sites.summary_gene.tsv", project=config['project'], group=samples["group"].dropna().unique()),
# # expand("{project}/groups/{group}/sig_xlsites/{group}.group.sig_sites.summary_type.tsv", project=config['project'], group=samples["group"].dropna().unique()),
# # expand("{project}/groups/{group}/sig_xlsites/{group}.group.sig_sites.summary_subtype.tsv", project=config['project'], group=samples["group"].dropna().unique()),
#
# expand("{project}/groups/{group}/clusters/{group}.group.clusters.bed", project=config['project'],
#        group=samples["group"].dropna().unique()),


#-------------------
# expand("{project}/groups/{group}/xlsites/{group}.group.unique.xl.bed", project=config['project'], group=samples["group"].dropna().unique()),
# expand("{project}/groups/{group}/xlsites/{group}.group.unique.xl.annotated_group_biotype.tab", project=config['project'], group=samples["group"].dropna().unique()),
# expand("{project}/groups/{group}/xlsites/{group}.group.unique.xl.annotated_group_gene_id.tab", project=config['project'], group=samples["group"].dropna().unique()),
# expand("{project}/groups/{group}/xlsites/{group}.group.unique.xl.summary_gene.tsv", project=config['project'], group=samples["group"].dropna().unique()),
# # expand("{project}/groups/{group}/xlsites/{group}.group.unique.xl.summary_type.tsv", project=config['project'], group=samples["group"].dropna().unique()),
# # expand("{project}/groups/{group}/xlsites/{group}.group.unique.xl.summary_subtype.tsv", project=config['project'], group=samples["group"].dropna().unique()),
# expand("{project}/groups/{group}/xlsites/{group}.group.unique.xl.bedgraph", project=config['project'], group=samples["group"].dropna().unique()),
# expand("{project}/groups/{group}/rnamaps/", project=config['project'], group=samples["group"].dropna().unique()),

# expand("{project}/groups/{group}/sig_xlsites/{group}.group.sig_sites.bed", project=config['project'], group=samples["group"].dropna().unique()),
# expand("{project}/groups/{group}/sig_xlsites/{group}.group.sig_sites.annotated.biotype.tab", project=config['project'], group=samples["group"].dropna().unique()),
# expand("{project}/groups/{group}/sig_xlsites/{group}.group.sig_sites.annotated.gene_id.tab", project=config['project'], group=samples["group"].dropna().unique()),
# expand("{project}/groups/{group}/sig_xlsites/{group}.group.sig_sites.summary_gene.tsv", project=config['project'], group=samples["group"].dropna().unique()),
# expand("{project}/groups/{group}/sig_xlsites/{group}.group.sig_sites.summary_type.tsv", project=config['project'], group=samples["group"].dropna().unique()),
# expand("{project}/groups/{group}/sig_xlsites/{group}.group.sig_sites.summary_subtype.tsv", project=config['project'], group=samples["group"].dropna().unique()),

#==============================================================================#
#                       SHA Check
#==============================================================================#


# print ("ALL_DEMULTIPLEXED", expand("demultiplexed/demux_{barcode}.fastq.gz", barcode=samples.index))
#  demultiplexed/demux_NNNNGTAACNNN_NNATT.fastq.gz demultiplexed/demux_NNNNGTAACNNN_NNAGG.fastq.gz demultiplexed/demux_NNNNGTAACNNN_NNTTA.fastq.gz demultiplexed/demux_NNNNGTAACNNN_NNTGC.fastq.gz demultiplexed/demux_NNNNGTAACNNN_NNCTG.fastq.gz demultiplexed/demux_NNNNGTAACNNN_NNCGT.fastq.gz demultiplexed/demux_NNNNGTAACNNN_NNGTC.fastq.gz demultiplexed/demux_NNNNGTAACNNN_NNGGA.fastq.gz demultiplexed/demux_NNNNCCGGANNN.fastq.gz demultiplexed/demux_NNNCTGCNN.fastq.gz
# cat  demultiplexed/demux_NNNNGTAACNNN_NNATT.fastq.gz demultiplexed/demux_NNNNGTAACNNN_NNAGG.fastq.gz demultiplexed/demux_NNNNGTAACNNN_NNTTA.fastq.gz demultiplexed/demux_NNNNGTAACNNN_NNTGC.fastq.gz demultiplexed/demux_NNNNGTAACNNN_NNCTG.fastq.gz demultiplexed/demux_NNNNGTAACNNN_NNCGT.fastq.gz demultiplexed/demux_NNNNGTAACNNN_NNGTC.fastq.gz demultiplexed/demux_NNNNGTAACNNN_NNGGA.fastq.gz demultiplexed/demux_NNNNCCGGANNN.fastq.gz demultiplexed/demux_NNNCTGCNN.fastq.gz | md5
# 3ef0e4c8a7628d7333df3cc315f498ce
# shasum demultiplexed/demux_NNNNGTAACNNN_NNATT.fastq.gz demultiplexed/demux_NNNNGTAACNNN_NNAGG.fastq.gz demultiplexed/demux_NNNNGTAACNNN_NNTTA.fastq.gz demultiplexed/demux_NNNNGTAACNNN_NNTGC.fastq.gz demultiplexed/demux_NNNNGTAACNNN_NNCTG.fastq.gz demultiplexed/demux_NNNNGTAACNNN_NNCGT.fastq.gz demultiplexed/demux_NNNNGTAACNNN_NNGTC.fastq.gz demultiplexed/demux_NNNNGTAACNNN_NNGGA.fastq.gz demultiplexed/demux_NNNNCCGGANNN.fastq.gz demultiplexed/demux_NNNCTGCNN.fastq.gz > demultiplex_md5.txt
# shasum -c demultiplex_md5.txt

# "qc/fastqc/raw_fastq_file_fastqc.html"
#

# "{project}/qc/fastqc/{barcode}_trimmed_fastqc.html"


#xlsites_output=list(xlsites_output)
rule shasum_create:
    input:
        demultiplex=expand("demultiplexed/demux_{barcode}.fastq.gz", barcode=samples.index),
        qc=["qc/fastqc/raw_fastq_file_fastqc.html",
            expand("{project}/qc/fastqc/{barcode}_fastqc.html", project=config['project'], barcode=samples.index),
            expand("{project}/qc/fastqc/{barcode}_trimmed_fastqc.html", project=config['project'], barcode=samples.index)],
        genome=[expand("{genomes_path}/{genome}/{genome}.fa.gz", genome=samples["mapto"].unique(), genomes_path=config['genomes_path']),
                expand("{genomes_path}/{genome}/{genome}.fa.gz.fai", genome=samples["mapto"].unique(), genomes_path=config['genomes_path']),
                expand("{genomes_path}/{genome}/{genome}.gtf.gz", genome=samples["mapto"].unique(), genomes_path=config['genomes_path']),
                expand("{genomes_path}/{genome}/segment/{genome}_segment.gtf", genome=samples["mapto"].unique(), genomes_path=config['genomes_path']),
                expand("{genomes_path}/{genome}/segment/landmarks.bed.gz", genome=samples["mapto"].unique(), genomes_path=config['genomes_path'])],
        mapped_reads=expand("{project}/mapped/{barcode}/Aligned.sortedByCoord.out.bam", project=config['project'], barcode=samples.index),
        xlsites=all_xlsites,
        group=all_group,

    output:
        demultiplex_shasum_file="{project}/shasum_files/demultiplex_shasum_file.txt",
        qc_shasum_file="{project}/shasum_files/qc_shasum_file.txt",
        genome_shasum_file="{project}/shasum_files/genome_shasum_file.txt",
        mapped_reads_shasum_file="{project}/shasum_files/mapped_reads_shasum_file.txt",
        xlsites_shasum_file="{project}/shasum_files/xlsites_shasum_file.txt",
        group_shasum_file="{project}/shasum_files/group_shasum_file.txt",
    shell:
        """
        shasum {input.demultiplex} > {output.demultiplex_shasum_file}
        shasum {input.qc} > {output.qc_shasum_file}
        shasum {input.genome} > {output.genome_shasum_file}
        shasum {input.mapped_reads} > {output.mapped_reads_shasum_file}
        shasum {input.xlsites} > {output.xlsites_shasum_file}
        shasum {input.group} > {output.group_shasum_file}
        """


rule shasum_check:
    input:
        demultiplex_shasum_file="{project}/shasum_files/demultiplex_shasum_file.txt",
        qc_shasum_file="{project}/shasum_files/qc_shasum_file.txt",
        genome_shasum_file="{project}/shasum_files/genome_shasum_file.txt",
        mapped_reads_shasum_file="{project}/shasum_files/mapped_reads_shasum_file.txt",
        xlsites_shasum_file="{project}/shasum_files/xlsites_shasum_file.txt",
        group_shasum_file="{project}/shasum_files/group_shasum_file.txt",
    output:
        demultiplex_shasum_check = "{project}/shasum_files/demultiplex_shasum_check.txt",
        qc_shasum_check="{project}/shasum_files/qc_shasum_check.txt",
        genome_shasum_check="{project}/shasum_files/genome_shasum_check.txt",
        mapped_reads_shasum_check="{project}/shasum_files/mapped_reads_shasum_check.txt",
        xlsites_shasum_check="{project}/shasum_files/xlsites_shasum_check.txt",
        group_shasum_check="{project}/shasum_files/group_shasum_check.txt",
    shell:
        """
        echo "Checking test files integrity and reproducibility."
        echo "--------------------------------------------------"
        shasum --check {input.demultiplex_shasum_file} 2>&1 | tee {output.demultiplex_shasum_check}
        shasum --check {input.qc_shasum_file} 2>&1 | tee {output.qc_shasum_check}
        shasum --check {input.genome_shasum_file} 2>&1 | tee {output.genome_shasum_check}
        shasum --check {input.mapped_reads_shasum_file} 2>&1 | tee {output.mapped_reads_shasum_check}
        shasum --check {input.xlsites_shasum_file} 2>&1 | tee {output.xlsites_shasum_check}
        shasum --check {input.group_shasum_file} 2>&1 | tee {output.group_shasum_check}
        """




#==============================================================================#
#                      Kmers
#==============================================================================#



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