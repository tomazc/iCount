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



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~* Final outputs *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

##### load rules #####
include: "rules/common.smk"
include: "rules/demultiplex.smk"
include: "rules/qc.smk"
include: "rules/prepare_genome.smk"
include: "rules/map_reads.smk"
include: "rules/bedgraph_UCSC.smk"



def all_input(wildcards):
    """
    Function defining all requested inputs for the rule all.
    """
    final_output = []

    if config["bedgraph_UCSC"]:
        final_output.extend(
            expand("{project}/xlsites/{barcode}/{barcode}.unique.xl.UCSC.bedgraph", project=config['project'], barcode=samples.index)
        )

    # if config["completeness_output"] == "minimum":
    #     final_output.extend(
    #         expand("{project}/xlsites/{barcode}/{barcode}.unique.xl.UCSC.bedgraph", project=config['project'], barcode=samples.index)
    #     )
    if config["create_integrity_test"]:
        final_output.extend(expand("{project}/shasum_files/demultiplex_shasum_file.txt", project=config['project']))
        final_output.extend(expand("{project}/shasum_files/qc_shasum_file.txt", project=config['project']))
        final_output.extend(expand("{project}/shasum_files/genome_shasum_file.txt", project=config['project']))

    if config["integrity_test_check"]:
        final_output.extend(expand("{project}/shasum_files/demultiplex_shasum_check.txt", project=config['project']))
        final_output.extend(expand("{project}/shasum_files/qc_shasum_check.txt", project=config['project']))
        final_output.extend(expand("{project}/shasum_files/genome_shasum_check.txt", project=config['project']))

    return final_output


##### target rules #####

localrules: all

rule all:
    input:
        "demultiplexed/demux_nomatch5.fastq.gz",

        expand("{project}/mapped/{barcode}/Aligned.sortedByCoord.out.bam", project=config['project'], barcode=samples.index),

        expand("{project}/xlsites/{barcode}/{barcode}.unique.xl.bed", project=config['project'], barcode=samples.index),
        expand("{project}/xlsites/{barcode}/{barcode}.unique.xl.summary_gene.tsv", project=config['project'], barcode=samples.index),
        expand("{project}/xlsites/{barcode}/{barcode}.unique.xl.annotated_sites_biotype.tab", project=config['project'], barcode=samples.index),
        expand("{project}/xlsites/{barcode}/{barcode}.unique.xl.annotated_sites_gene_id.tab", project=config['project'], barcode=samples.index),
        expand("{project}/xlsites/{barcode}/{barcode}.unique.xl.bedgraph", project=config['project'], barcode=samples.index),
        #expand("{project}/xlsites/{barcode}/{barcode}.unique.xl.UCSC.bedgraph", project=config['project'], barcode=samples.index),
        expand("{project}/xlsites/{barcode}/rnamaps/", project=config['project'], barcode=samples.index),

        expand("{project}/sig_xlsites/{barcode}/{barcode}.sig_sites.bed", project=config['project'], barcode=samples.index),
        expand("{project}/sig_xlsites/{barcode}/{barcode}.sig_sites.summary_gene.tsv", project=config['project'], barcode=samples.index),
        expand("{project}/sig_xlsites/{barcode}/{barcode}.sig_sites.annotated.biotype.tab", project=config['project'], barcode=samples.index),
        expand("{project}/sig_xlsites/{barcode}/{barcode}.sig_sites.annotated.gene_id.tab", project=config['project'], barcode=samples.index),

        expand("{project}/clusters/{barcode}/{barcode}.clusters.bed", project=config['project'], barcode=samples.index),

        expand("{project}/groups/{group}/xlsites/{group}.group.unique.xl.bed", project=config['project'], group=samples["group"].dropna().unique()),
        expand("{project}/groups/{group}/xlsites/{group}.group.unique.xl.annotated_group_biotype.tab", project=config['project'], group=samples["group"].dropna().unique()),
        expand("{project}/groups/{group}/xlsites/{group}.group.unique.xl.annotated_group_gene_id.tab", project=config['project'], group=samples["group"].dropna().unique()),
        expand("{project}/groups/{group}/xlsites/{group}.group.unique.xl.summary_gene.tsv", project=config['project'], group=samples["group"].dropna().unique()),
        # expand("{project}/groups/{group}/xlsites/{group}.group.unique.xl.summary_type.tsv", project=config['project'], group=samples["group"].dropna().unique()),
        # expand("{project}/groups/{group}/xlsites/{group}.group.unique.xl.summary_subtype.tsv", project=config['project'], group=samples["group"].dropna().unique()),
        expand("{project}/groups/{group}/xlsites/{group}.group.unique.xl.bedgraph", project=config['project'], group=samples["group"].dropna().unique()),
        expand("{project}/groups/{group}/rnamaps/", project=config['project'], group=samples["group"].dropna().unique()),

        # expand("{project}/groups/{group}/sig_xlsites/{group}.group.sig_sites.bed", project=config['project'], group=samples["group"].dropna().unique()),
        # expand("{project}/groups/{group}/sig_xlsites/{group}.group.sig_sites.annotated.biotype.tab", project=config['project'], group=samples["group"].dropna().unique()),
        # expand("{project}/groups/{group}/sig_xlsites/{group}.group.sig_sites.annotated.gene_id.tab", project=config['project'], group=samples["group"].dropna().unique()),
        # expand("{project}/groups/{group}/sig_xlsites/{group}.group.sig_sites.summary_gene.tsv", project=config['project'], group=samples["group"].dropna().unique()),
        # expand("{project}/groups/{group}/sig_xlsites/{group}.group.sig_sites.summary_type.tsv", project=config['project'], group=samples["group"].dropna().unique()),
        # expand("{project}/groups/{group}/sig_xlsites/{group}.group.sig_sites.summary_subtype.tsv", project=config['project'], group=samples["group"].dropna().unique()),

        expand("{project}/groups/{group}/clusters/{group}.group.clusters.bed", project=config['project'], group=samples["group"].dropna().unique()),
        all_input

# old rule all not used - remove
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
    output:
        demultiplex_shasum_file="{project}/shasum_files/demultiplex_shasum_file.txt",
        qc_shasum_file="{project}/shasum_files/qc_shasum_file.txt",
        genome_shasum_file="{project}/shasum_files/genome_shasum_file.txt",
        mapped_reads_shasum_file="{project}/shasum_files/mapped_reads_shasum_file.txt",
    shell:
        """
        shasum {input.demultiplex} > {output.demultiplex_shasum_file}
        shasum {input.qc} > {output.qc_shasum_file}
        shasum {input.genome} > {output.genome_shasum_file}
        shasum {input.mapped_reads} > {output.mapped_reads_shasum_file}
        """

rule shasum_check:
    input:
        demultiplex_shasum_file="{project}/shasum_files/demultiplex_shasum_file.txt",
        qc_shasum_file="{project}/shasum_files/qc_shasum_file.txt",
        genome_shasum_file="{project}/shasum_files/genome_shasum_file.txt",
        mapped_reads_shasum_file="{project}/shasum_files/mapped_reads_shasum_file.txt",
    output:
        demultiplex_shasum_check = "{project}/shasum_files/demultiplex_shasum_check.txt",
        qc_shasum_check="{project}/shasum_files/qc_shasum_check.txt",
        genome_shasum_check="{project}/shasum_files/genome_shasum_check.txt",
        mapped_reads_shasum_check="{project}/shasum_files/mapped_reads_shasum_check.txt",
    shell:
        """
        echo "Checking test files integrity and reproducibility."
        echo "--------------------------------------------------"
        shasum --check {input.demultiplex_shasum_file} 2>&1 | tee {output.demultiplex_shasum_check}
        shasum --check {input.qc_shasum_file} 2>&1 | tee {output.qc_shasum_check}
        shasum --check {input.genome_shasum_file} 2>&1 | tee {output.genome_shasum_check}
        shasum --check {input.mapped_reads_shasum_file} 2>&1 | tee {output.mapped_reads_shasum_check}
        """




#==============================================================================#
#     Create crosslink sites, annotate and obtain summaries
#==============================================================================#


rule xlsites:
    input:
        "{project}/mapped/{barcode}/Aligned.sortedByCoord.out.bam"
    output:
        unique_bed="{project}/xlsites/{barcode}/{barcode}.unique.xl.bed",
        multimapped_bed="{project}/xlsites/{barcode}/{barcode}.multimapped.xl.bed",
        skipped_bam="{project}/xlsites/{barcode}/{barcode}.skipped.xl.bam",
    benchmark:
        "{project}/benchmarks/{barcode}.xlsites.benchmark.tab"
        # repeat("benchmarks/{barcode}.xlsites.benchmark.tab", 3)
    params:
        group_by=config['group_by'],
        quant = config['quant'],
        mismatches = config['mismatches'],
        mapq_th = config['mapq_th'],
        multimax = config['multimax'],
        gap_th = config['gap_th'],
        ratio_th = config['ratio_th'],
        max_barcodes = config['max_barcodes'],
    log:
        "{project}/logs/xlsites/{barcode}.xlsites.log"
    shell:
        """
        iCount xlsites --group_by {params.group_by} --quant {params.quant} -mis {params.mismatches} --mapq_th {params.mapq_th} \
        --multimax {params.multimax} --gap_th {params.gap_th} --ratio_th {params.ratio_th} --max_barcodes {params.max_barcodes} \
        {input} {output.unique_bed} {output.multimapped_bed} {output.skipped_bam} -M {log}

        """

def bedgraph_description(wildcards):
    return ("{project}_{sample_name}_{protein}_{method}_{mapto}".format(project=config['project'], sample_name=samples.loc[wildcards.barcode, "sample_name"], mapto=samples.loc[wildcards.barcode, "mapto"], method=samples.loc[wildcards.barcode, "method"],	protein=samples.loc[wildcards.barcode, "protein"],	cells_tissue=samples.loc[wildcards.barcode, "cells_tissue"],	condition=samples.loc[wildcards.barcode, "condition"],))



rule bedgraph:
    input:
        xlsites="{project}/xlsites/{barcode}/{barcode}.unique.xl.bed",
    output:
        bedgraph="{project}/xlsites/{barcode}/{barcode}.unique.xl.bedgraph",
    params:
        name=bedgraph_description,
        description=bedgraph_description,
        visibility="full",
        priority="20",
        color="120,101,172",
        alt_color="200,120,59",
        max_height_pixels="100:50:0",
    run:
        shell("iCount bedgraph --name \"{params.name}.unique.xl.bedgraph\" --description \"{params.description}\" --visibility \"{params.visibility}\" --priority \"{params.priority}\" --color \"{params.color}\" --alt_color \"{params.alt_color}\" --max_height_pixels \"{params.max_height_pixels}\" {input.xlsites} {output.bedgraph}")


rule annotate_xlsites:
    input:
        xlsites = "{project}/xlsites/{barcode}/{barcode}.unique.xl.bed"
    output:
        biotype="{project}/xlsites/{barcode}/{barcode}.unique.xl.annotated_sites_biotype.tab",
        gene_id="{project}/xlsites/{barcode}/{barcode}.unique.xl.annotated_sites_gene_id.tab",
    params:
        templates_dir = get_templates_dir,
        segment = get_segment_path,
        #out_dir = "{project}/annotated/",
    shell:
        """
        iCount annotate --subtype biotype {params.segment} {input.xlsites} {output.biotype}
        iCount annotate --subtype gene_id {params.segment} {input.xlsites} {output.gene_id}
        """

rule summary:
    input:
        xlsites="{project}/xlsites/{barcode}/{barcode}.unique.xl.bed"
    output:
        gene="{project}/xlsites/{barcode}/{barcode}.unique.xl.summary_gene.tsv",
        type="{project}/xlsites/{barcode}/{barcode}.unique.xl.summary_type.tsv",
        subtype="{project}/xlsites/{barcode}/{barcode}.unique.xl.summary_subtype.tsv"
    params:
        templates_dir=get_templates_dir,
        segment = get_segment_path,
        out_dir="{project}/xlsites/{barcode}/",
        rename_gene="{project}/xlsites/{barcode}/summary_gene.tsv",
        rename_type="{project}/xlsites/{barcode}/summary_type.tsv",
        rename_subtype="{project}/xlsites/{barcode}/summary_subtype.tsv",
    shell:
        """
        iCount summary --templates_dir {params.templates_dir} {params.segment} {input.xlsites} {params.out_dir}
        mv {params.rename_gene} {output.gene}
        mv {params.rename_type} {output.type}
        mv {params.rename_subtype} {output.subtype}
        """


def get_xlsites_landmark_path(wildcards):
    return ("{0}/{1}/segment/landmarks.bed.gz".format(config['genomes_path'], samples.loc[wildcards.barcode, "mapto"]))


rule RNAmaps:
    input:
        xlsites="{project}/xlsites/{barcode}/{barcode}.unique.xl.bed",
        landmarks=get_xlsites_landmark_path,
    output:
        directory("{project}/xlsites/{barcode}/rnamaps/")
    params:
        plot_type=config['plot_type'],
        top_n=config['top_n'],
        smoothing=config['smoothing'],
        nbins=config['nbins'],
        binsize=config['binsize'],
        colormap=config['colormap'],
        imgfmt=config['imgfmt'],
    shell:
        """
        iCount rnamaps {input.xlsites} {input.landmarks} --outdir {output} --top_n {params.top_n} \
        --smoothing {params.smoothing} --colormap {params.colormap} --imgfmt {params.imgfmt} --nbins {params.nbins} 
        # --binsize {params.binsize}
        """

#==============================================================================#
#       Create significant crosslink sites, annotate and obtain summaries
#==============================================================================#

rule sig_xlsites:
    input:
        xlsites="{project}/xlsites/{barcode}/{barcode}.unique.xl.bed",
        segment_file=get_segment_path
    output:
        sigxls="{project}/sig_xlsites/{barcode}/{barcode}.sig_sites.bed",
        scores="{project}/sig_xlsites/{barcode}/{barcode}.scores.tsv"
    benchmark:
        "{project}/benchmarks/{barcode}.sig_xlsites.benchmark.txt"
    shell:
        """
        iCount peaks {input.segment_file} {input.xlsites} {output.sigxls} --scores {output.scores}
        """


def is_empty(fname):
    print(fname + " file is empty.")
    return os.stat(str(fname)).st_size == 0


# Create a new empty file.
def createNewFile(fname):
    file_object = open(fname, 'w')
    # file_object.write('File is created.')
    print(fname + " has been created. ")


rule annotate_sig_xlsites:
    input:
        sig_xlsites = "{project}/sig_xlsites/{barcode}/{barcode}.sig_sites.bed"
    output:
        biotype="{project}/sig_xlsites/{barcode}/{barcode}.sig_sites.annotated.biotype.tab",
        gene_id="{project}/sig_xlsites/{barcode}/{barcode}.sig_sites.annotated.gene_id.tab"
    params:
        templates_dir = get_templates_dir,
        segment = get_segment_path,
        #out_dir = "{project}/sig_xlsites/{barcode}/",
    run:
        if is_empty(input.sig_xlsites):
            print ("File", input.sig_xlsites, "is empty. Creating empty output files:", output.biotype, output.gene_id, \
                   " to continue snakemake pipeline")
            createNewFile(output.biotype)
            createNewFile(output.gene_id)
        else:
            shell("iCount annotate --subtype biotype {params.segment} {input.sig_xlsites} {output.biotype}")
            shell("iCount annotate --subtype gene_id {params.segment} {input.sig_xlsites} {output.gene_id}")


rule summary_sig:
    input:
        sig_xlsites = "{project}/sig_xlsites/{barcode}/{barcode}.sig_sites.bed"
    output:
        gene = "{project}/sig_xlsites/{barcode}/{barcode}.sig_sites.summary_gene.tsv",
        type = "{project}/sig_xlsites/{barcode}/{barcode}.sig_sites.summary_type.tsv",
        subtype = "{project}/sig_xlsites/{barcode}/{barcode}.sig_sites.summary_subtype.tsv",
    params:
        templates_dir = get_templates_dir,
        segment = get_segment_path,
        out_dir = "{project}/sig_xlsites/{barcode}/",
        rename_gene = "{project}/sig_xlsites/{barcode}/summary_gene.tsv",
        rename_type = "{project}/sig_xlsites/{barcode}/summary_type.tsv",
        rename_subtype = "{project}/sig_xlsites/{barcode}/summary_subtype.tsv",
    run:
        if is_empty(input.sig_xlsites):
            print ("File", input.sig_xlsites, "is empty. Creating empty output file:", output.gene,
                   " to continue snakemake pipeline")
            createNewFile(output.gene)
            createNewFile(output.type)
            createNewFile(output.subtype)
        else:
            shell("iCount summary --templates_dir {params.templates_dir} {params.segment} {input.sig_xlsites} {params.out_dir}")
            shell("mv {params.rename_gene} {output.gene}")
            shell("mv {params.rename_type} {output.type}")
            shell("mv {params.rename_subtype} {output.subtype}")


#==============================================================================#
#             Create clusters
#==============================================================================#

rule clusters:
    input:
        xlsites="{project}/xlsites/{barcode}/{barcode}.unique.xl.bed",
        sig_xlsites="{project}/sig_xlsites/{barcode}/{barcode}.sig_sites.bed",
    output:
        clusters="{project}/clusters/{barcode}/{barcode}.clusters.bed",
    benchmark:
        "{project}/benchmarks/{barcode}.clusters.benchmark.txt"
    params:
        distance=config['distance'],
        slop=config['slop'],
    run:
        if is_empty(input.xlsites) or is_empty(input.sig_xlsites):
            print ("File", input.sig_xlsites, "is empty. Creating empty output file:", output.clusters, \
                   " to continue snakemake pipeline")
            createNewFile(output.clusters)
        else:
            shell("iCount clusters --dist {params.distance} --slop {params.slop} {input.xlsites} {input.sig_xlsites} \
            {output.clusters}")




#==============================================================================#
#                       Group analysis
#==============================================================================#

def files2group(wildcards):
    barcode_group = samples.loc[samples['group'] == wildcards.group, "full_barcode"].values[0:].tolist()
    xlsites_list=[]
    for barcode in barcode_group:
        xlsites_list.append("{project}/xlsites/{barcode}/{barcode}.unique.xl.bed".format(project=config['project'], barcode=barcode))
    return (xlsites_list)


rule group:
    input:
        files2group,
    output:
        group="{project}/groups/{group}/xlsites/{group}.group.unique.xl.bed",
    benchmark:
        "{project}/benchmarks/{group}.group.unique.xl.benchmark.txt"
    shell:
        """
        iCount group {output.group} {input}
        """

def bedgraph_group_description(wildcards):
    return ("{project}_group_{group}_{protein}_{method}_{mapto}".format(project=config['project'], group=wildcards.group, mapto=samples.loc[samples['group'] == wildcards.group, "mapto"].unique()[0], method=samples.loc[samples['group'] == wildcards.group, "method"].unique()[0],	protein=samples.loc[samples['group'] == wildcards.group, "protein"].unique()[0],	cells_tissue=samples.loc[samples['group'] == wildcards.group, "cells_tissue"].unique()[0],	condition=samples.loc[samples['group'] == wildcards.group, "condition"].unique()[0]))



rule group_bedgraph:
    input:
        group_xlsites="{project}/groups/{group}/xlsites/{group}.group.unique.xl.bed",
    output:
        group_bedgraph="{project}/groups/{group}/xlsites/{group}.group.unique.xl.bedgraph",
    params:
        name=bedgraph_group_description,
        description=bedgraph_group_description,
        visibility="full",
        priority="20",
        color="120,101,172",
        alt_color="200,120,59",
        max_height_pixels="100:50:0",
    run:
        shell("iCount bedgraph --name \"{params.name}.group.unique.xl.bedgraph\" --description \"{params.description}\" --visibility \"{params.visibility}\" --priority \"{params.priority}\" --color \"{params.color}\" --alt_color \"{params.alt_color}\" --max_height_pixels \"{params.max_height_pixels}\" {input.group_xlsites} {output.group_bedgraph}")



def get_group_segment_path(wildcards):
    return ("{0}/{1}/segment/{1}_segment.gtf".format(config['genomes_path'], samples.loc[samples['group'] == wildcards.group, "mapto"].unique()[0]))

def get_group_templates_dir(wildcards):
    return ("{0}/{1}/segment/".format(config['genomes_path'], samples.loc[samples['group'] == wildcards.group, "mapto"].unique()[0]))

rule annotate_group_xlsites:
    input:
        group_xlsites="{project}/groups/{group}/xlsites/{group}.group.unique.xl.bed"
    output:
        group_biotype="{project}/groups/{group}/xlsites/{group}.group.unique.xl.annotated_group_biotype.tab",
        group_gene_id="{project}/groups/{group}/xlsites/{group}.group.unique.xl.annotated_group_gene_id.tab",
    params:
        templates_dir = get_group_templates_dir,
        segment = get_group_segment_path,
        #out_dir = "{project}/annotated/",
    shell:
        """
        iCount annotate --subtype biotype {params.segment} {input.group_xlsites} {output.group_biotype}
        iCount annotate --subtype gene_id {params.segment} {input.group_xlsites} {output.group_gene_id}
        """

rule summary_group:
    input:
        group_xlsites="{project}/groups/{group}/xlsites/{group}.group.unique.xl.bed"
    output:
        group_gene="{project}/groups/{group}/xlsites/{group}.group.unique.xl.summary_gene.tsv",
        group_type="{project}/groups/{group}/xlsites/{group}.group.unique.xl.summary_type.tsv",
        group_subtype="{project}/groups/{group}/xlsites/{group}.group.unique.xl.summary_subtype.tsv"
    params:
        templates_dir=get_group_templates_dir,
        segment = get_group_segment_path,
        out_dir="{project}/groups/{group}/xlsites/",
        group_rename_gene="{project}/groups/{group}/xlsites/summary_gene.tsv",
        group_rename_type="{project}/groups/{group}/xlsites/summary_type.tsv",
        group_rename_subtype="{project}/groups/{group}/xlsites/summary_subtype.tsv",
    shell:
        """
        iCount summary --templates_dir {params.templates_dir} {params.segment} {input.group_xlsites} {params.out_dir}
        mv {params.group_rename_gene} {output.group_gene}
        mv {params.group_rename_type} {output.group_type}
        mv {params.group_rename_subtype} {output.group_subtype}
        """

def get_group_segment_path(wildcards):
    return ("{0}/{1}/segment/{1}_segment.gtf".format(config['genomes_path'], samples.loc[samples['group'] == wildcards.group, "mapto"].unique()[0]))

rule group_sig_xlsites:
    input:
        group_xlsites = "{project}/groups/{group}/xlsites/{group}.group.unique.xl.bed",
        segment_file=get_group_segment_path,
    output:
        group_sigxls="{project}/groups/{group}/sig_xlsites/{group}.group.sig_sites.bed",
        group_scores="{project}/groups/{group}/sig_xlsites/{group}.group.scores.tsv"
    benchmark:
        "{project}/benchmarks/{group}.group.sig_xlsites.benchmark.txt"
    shell:
        """
        iCount peaks {input.segment_file} {input.group_xlsites} {output.group_sigxls} --scores {output.group_scores}
        """

rule annotate_group_sig_xlsites:
    input:
        group_sig_xlsites="{project}/groups/{group}/sig_xlsites/{group}.group.sig_sites.bed"
    output:
        group_biotype="{project}/groups/{group}/sig_xlsites/{group}.group.sig_sites.annotated.biotype.tab",
        group_gene_id="{project}/groups/{group}/sig_xlsites/{group}.group.sig_sites.annotated.gene_id.tab"
    params:
        templates_dir = get_group_templates_dir,
        segment = get_group_segment_path,
        #out_dir = "{project}/sig_xlsites/{barcode}/",
    run:
        if is_empty(input.group_sig_xlsites):
            print ("File", input.group_sig_xlsites, "is empty. Creating empty output files:", output.group_biotype, \
                   output.group_gene_id, " to continue snakemake pipeline")
            createNewFile(output.group_biotype)
            createNewFile(output.group_gene_id)
        else:
            shell("iCount annotate --subtype biotype {params.segment} {input.group_sig_xlsites} {output.group_biotype}")
            shell("iCount annotate --subtype gene_id {params.segment} {input.group_sig_xlsites} {output.group_gene_id}")


rule summary_group_sig:
    input:
        group_sig_xlsites="{project}/groups/{group}/sig_xlsites/{group}.group.sig_sites.bed"
    output:
        group_gene="{project}/groups/{group}/sig_xlsites/{group}.group.sig_sites.summary_gene.tsv",
        group_type="{project}/groups/{group}/sig_xlsites/{group}.group.sig_sites.summary_type.tsv",
        group_subtype="{project}/groups/{group}/sig_xlsites/{group}.group.sig_sites.summary_subtype.tsv",
    params:
        templates_dir = get_group_templates_dir,
        segment = get_group_segment_path,
        out_dir = "{project}/groups/{group}/sig_xlsites/",
        group_rename_gene = "{project}/groups/{group}/sig_xlsites/summary_gene.tsv",
        group_rename_type = "{project}/groups/{group}/sig_xlsites/summary_type.tsv",
        group_rename_subtype = "{project}/groups/{group}/sig_xlsites/summary_subtype.tsv",
    run:
        if is_empty(input.group_sig_xlsites):
            print ("File", input.group_sig_xlsites, "is empty. Creating empty output file:", output.group_gene, \
                   output.group_type , output.group_subtype, " to continue snakemake pipeline")
            createNewFile(output.group_gene)
            createNewFile(output.group_type)
            createNewFile(output.group_subtype)
        else:
            shell("iCount summary --templates_dir {params.templates_dir} {params.segment} {input.sig_xlsites} {params.out_dir}")
            shell("mv {params.group_rename_gene} {output.group_gene}")
            shell("mv {params.group_rename_type} {output.group_type}")
            shell("mv {params.group_rename_subtype} {output.group_subtype}")


def get_group_landmark_path(wildcards):
    return ("{0}/{1}/segment/landmarks.bed.gz".format(config['genomes_path'], samples.loc[samples['group'] == wildcards.group, "mapto"].unique()[0]))

rule group_clusters:
    input:
        group_xlsites="{project}/groups/{group}/xlsites/{group}.group.unique.xl.bed",
        group_sigxls="{project}/groups/{group}/sig_xlsites/{group}.group.sig_sites.bed",
    output:
        group_clusters="{project}/groups/{group}/clusters/{group}.group.clusters.bed",
    benchmark:
        "{project}/benchmarks/{group}.group.clusters.benchmark.txt"
    params:
        distance=config['distance'],
        slop=config['slop'],
    run:
        if is_empty(input.group_xlsites) or is_empty(input.group_sigxls):
            print ("File", input.group_sigxls, "is empty. Creating empty output file:", output.group_clusters, \
                   " to continue snakemake pipeline")
            createNewFile(output.group_clusters)
        else:
            shell("iCount clusters --dist {params.distance} --slop {params.slop} {input.group_xlsites} \
            {input.group_sigxls} {output.group_clusters}")


rule RNAmaps_group:
    input:
        group_xlsites="{project}/groups/{group}/xlsites/{group}.group.unique.xl.bed",
        landmarks=get_group_landmark_path,
    output:
        directory("{project}/groups/{group}/rnamaps/")
    params:
        plot_type=config['plot_type'],
        top_n=config['top_n'],
        smoothing=config['smoothing'],
        nbins=config['nbins'],
        binsize=config['binsize'],
        colormap=config['colormap'],
        imgfmt=config['imgfmt'],
    shell:
        """
        iCount rnamaps {input.group_xlsites} {input.landmarks} --outdir {output} --top_n {params.top_n} --smoothing \
        {params.smoothing} --colormap {params.colormap} --imgfmt {params.imgfmt} --nbins {params.nbins} 
        # --binsize {params.binsize}
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