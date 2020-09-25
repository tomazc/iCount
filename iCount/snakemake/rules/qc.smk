
#==============================================================================#
#                       Read quality trimming and QC
#==============================================================================#

rule fastqc_raw:
    input:
        config['raw_fastq_file']
    output:
        html="qc/fastqc/raw_fastq_file_fastqc.html",
        zip="qc/fastqc/raw_fastq_file_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    wrapper:
        "0.38.0/bio/fastqc"


rule fastqc:
    input:
        "demultiplexed/demux_{barcode}.fastq.gz"
    output:
        html="{project}/qc/fastqc/{barcode}_fastqc.html",
        zip="{project}/qc/fastqc/{barcode}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    log:
        "{project}/logs/fastqc/{barcode}_fastqc.txt"
    wrapper:
        "0.36.0/bio/fastqc"


# Trimm reads
rule quality_trim:
    input:
        "demultiplexed/demux_{barcode}.fastq.gz"
    output:
        trimmed_reads="{project}/trimmed/demux_{barcode}_trimmed.fastq.gz",
        metrics="{project}/metrics/{barcode}_trimmed.txt"
    params:
        qual_trim=config['qual_trim'],
        minimum_length=config['minimum_length'],
        adapter=samples["adapter_3"].unique().tolist(),
        overlap=config['overlap'],
        untrimmed_output=config['untrimmed_output'],
        error_rate=config['error_rate'],
    log:
        "{project}/logs/trimmed/{barcode}_trimmed.txt"
    shell:
        """
        iCount cutadapt --qual_trim {params.qual_trim} --untrimmed_output {params.untrimmed_output} --minimum_length {params.minimum_length} --file_log 2 --file_logpath {log} --results_file {output.metrics} --reads_trimmed {output.trimmed_reads} {input} {params.adapter}
        """
        ## Unused parameters: --overlap {params.overlap}

rule fastqc_trimmed:
    input:
        "{project}/trimmed/demux_{barcode}_trimmed.fastq.gz"
    output:
        html="{project}/qc/fastqc/{barcode}_trimmed_fastqc.html",
        zip="{project}/qc/fastqc/{barcode}_trimmed_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    log:
        "{project}/logs/fastqc/{barcode}_trimmed_fastqc.log"
    wrapper:
        "0.36.0/bio/fastqc"

