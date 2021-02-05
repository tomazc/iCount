#==============================================================================#
#                       Demultiplex
#==============================================================================#
#### Check if one of the barcodes OR nomatch is not created/found

# def demux_out(wildcards):
#     demux_output = list()
#     #expand("demultiplexed/demux_{0}.fastq.gz".format(samples.loc[ "full_barcode"])))
#     demux_output.extend(expand("demultiplexed/demux_{barcode}.fastq.gz", barcode=samples.index))
#     demux_output.extend("demultiplexed/demux_nomatch5.fastq.gz")
#     print ("demux_output", demux_output)
#     return (demux_output)
#
#     # expand("demultiplexed/demux_{barcode}.fastq.gz", barcode=samples.index), \
#     # "demultiplexed/demux_nomatch5.fastq.gz"
#


rule demultiplex:
    input:
        fastq_file=config['raw_fastq_file'],
    output:
        expand("demultiplexed/demux_{barcode}.fastq.gz", barcode=samples.index),
        "demultiplexed/demux_nomatch5.fastq.gz"
    params:
        adapter3=samples["adapter_3"].unique().tolist(),
        all_5barcodes=samples["barcode_5"].tolist(),
        # all_3barcodes=samples["barcode_3"].tolist(),
        all_3barcodes = samples["barcode_3"].fillna(".").tolist(),
        dir=directory("demultiplexed"),
        barcode_mismatches=config['barcode_mismatches'],
        minimum_length=config['minimum_length'],
        min_adapter_overlap=config['min_adapter_overlap'],
    shell:
        """
        iCount demultiplex --mismatches {params.barcode_mismatches} --min_adapter_overlap {params.min_adapter_overlap} --minimum_length {params.minimum_length} {input.fastq_file} {params.adapter3} {params.all_5barcodes} --barcodes3 {params.all_3barcodes} --out_dir {params.dir}
        """






# rule move_demultiplex:
#     input:
#         directory("demultiplexed/")
#     output:
#         directory("{project}/demultiplexed/".format(project=config['project']))
#     run:
#         shutil.copytree(input, output)

# -M {log.metrics} 2> {log.log}
# log:
# metrics = "{project}/metrics/demultiplex_metrics.txt",
# log = "{project}/logs/demultiplex_metrics.txt"