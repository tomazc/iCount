#==============================================================================#
#                       Map reads
#==============================================================================#



rule map_reads:
    input:
        trimmed_reads="{project}/trimmed/demux_{barcode}_trimmed.fastq.gz",
        gtf = get_gtf_path,
    output:
        "{project}/mapped/{barcode}/Aligned.sortedByCoord.out.bam"
    params:
        star_index = directory(get_star_index_path),
        outdir=directory("{project}/mapped/{barcode}/"),
        multimax=config['multimax'],
    log:
        "{project}/logs/mapstar/{barcode}_mapstar.log"
    shell:
        """
        iCount mapstar --annotation {input.gtf} --multimax {params.multimax} \
        {input.trimmed_reads} {params.star_index} {params.outdir}
        """

