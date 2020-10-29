#==============================================================================#
#                       Map reads
#==============================================================================#



rule map_reads:
    input:
        trimmed_reads="{project}/trimmed/demux_{barcode}_qtrimmed.fastq.gz",
        gtf = get_gtf_path,
        star_index = get_star_index_folder,
    output:
        "{project}/mapped/{barcode}/Aligned.sortedByCoord.out.bam"
    params:
        outdir=directory("{project}/mapped/{barcode}/"),
        multimax=config['multimax'],
    log:
        "{project}/logs/mapstar/{barcode}_mapstar.log"
    shell:
        """
        iCount mapstar --annotation {input.gtf} --multimax {params.multimax} \
        {input.trimmed_reads} {input.star_index} {params.outdir}
        """

