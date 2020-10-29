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
        segment = get_segment_regions,
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
        segment = get_segment_regions,
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