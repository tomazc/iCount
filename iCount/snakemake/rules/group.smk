#==============================================================================#
#                       Group analysis
#==============================================================================#


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
            shell("iCount summary --templates_dir {params.templates_dir} {params.segment} {input.group_sig_xlsites} {params.out_dir}")
            shell("mv {params.group_rename_gene} {output.group_gene}")
            shell("mv {params.group_rename_type} {output.group_type}")
            shell("mv {params.group_rename_subtype} {output.group_subtype}")


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
