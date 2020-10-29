
#==============================================================================#
#                       Read quality trimming and QC
#==============================================================================#

# Check quality of raw reads
rule fastqc_raw:
    input:
        config['raw_fastq_file']
    output:
        html="{project}/qc/fastqc/raw_fastq_file_fastqc.html",
        zip="{project}/qc/fastqc/raw_fastq_file_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    wrapper:
        "0.38.0/bio/fastqc"

# Check quality of demultiplexed reads
rule fastqc:
    input:
        "demultiplexed/demux_{barcode}.fastq.gz"
    output:
        html="{project}/qc/fastqc/{barcode}_fastqc.html",
        zip="{project}/qc/fastqc/{barcode}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    log:
        "{project}/logs/fastqc/{barcode}_fastqc.txt"
    wrapper:
        "0.38.0/bio/fastqc"




# Trimm reads
rule quality_trim:
    input:
        "demultiplexed/demux_{barcode}.fastq.gz"
    output:
        qtrimmed_output = "{project}/trimmed/demux_{barcode}_qtrimmed.fastq.gz",
        metrics="{project}/metrics/{barcode}_qtrimmed.txt",
    params:
        trimmed_reads = "{project}/trimmed/demux_{barcode}_trimmed.fastq.gz",
        untrimmed_reads = "{project}/trimmed/demux_{barcode}_untrimmed.fastq.gz",
        qual_trim=config['qual_trim'],
        minimum_length=config['minimum_length'],
        adapter=samples["adapter_3"].unique().tolist(),
        overlap=config['overlap'],
        error_rate=config['error_rate'],
    log:
        "{project}/logs/trimmed/{barcode}_qtrimmed.txt"
    shell:
        """
        iCount cutadapt --qual_trim {params.qual_trim} --untrimmed_output {params.untrimmed_reads} --minimum_length {params.minimum_length} --file_log 2 --file_logpath {log} --results_file {output.metrics} --reads_trimmed {params.trimmed_reads} {input} {params.adapter}
        cat {params.trimmed_reads} {params.untrimmed_reads} > {output.qtrimmed_output}
        #rm {params.trimmed_reads}
        """
        ## Unused parameters: --overlap {params.overlap}

rule fastqc_trimmed:
    input:
        "{project}/trimmed/demux_{barcode}_qtrimmed.fastq.gz"
    output:
        html="{project}/qc/fastqc/{barcode}_qtrimmed_fastqc.html",
        zip="{project}/qc/fastqc/{barcode}_qtrimmed_trimmed_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    log:
        "{project}/logs/fastqc/{barcode}_qtrimmed_fastqc.log"
    wrapper:
        "0.38.0/bio/fastqc"

