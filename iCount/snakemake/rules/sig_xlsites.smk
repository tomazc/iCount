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
    params:
        features=config['features'],
        group_by=config['group_by_sig_xl'],
        merge_features=config['merge_features'],
        half_window=config['half_window'],
        fdr=config['fdr'],
        perms=config['perms'],
        rnd_seed=config['rnd_seed'],
        results_file="{project}/sig_xlsites/{barcode}/{barcode}.sig_sites_metrics.txt",
    log:
        file_logpath="{project}/logs/sig_xlsites/{barcode}.xlsites.log"
    benchmark:
        "{project}/benchmarks/{barcode}.sig_xlsites.benchmark.txt"
    shell:
        """
        iCount peaks {input.segment_file} {input.xlsites} {output.sigxls} --scores {output.scores} --features {params.features} --group_by {params.group_by} --half_window {params.half_window} --fdr {params.fdr} --perms {params.perms} --rnd_seed {params.rnd_seed} --file_logpath {log.file_logpath} --results_file {params.results_file} 
        """

# iCount peaks iCount_genomes/homo_sapiens/segment/regions.gtf.gz synthetic_chr20/xlsites/NNNNGTAACNNN_NNAGG/NNNNGTAACNNN_NNAGG.unique.xl.bed synthetic_chr20/sig_xlsites/NNNNGTAACNNN_NNAGG/NNNNGTAACNNN_NNAGG.sig_sites.bed --scores synthetic_chr20/sig_xlsites/NNNNGTAACNNN_NNAGG/NNNNGTAACNNN_NNAGG.scores.tsv --features gene --group_by gene_id --half_window 3 --fdr 0.05 --perms 100 --rnd_seed 42 --file_logpath synthetic_chr20/logs/sig_xlsites/NNNNGTAACNNN_NNAGG.xlsites.log --results_file synthetic_chr20/sig_xlsites/NNNNGTAACNNN_NNAGG/NNNNGTAACNNN_NNAGG.sig_sites_metrics.txt


rule annotate_sig_xlsites:
    input:
        sig_xlsites = "{project}/sig_xlsites/{barcode}/{barcode}.sig_sites.bed"
    output:
        biotype="{project}/sig_xlsites/{barcode}/{barcode}.sig_sites.annotated.biotype.tab",
        gene_id="{project}/sig_xlsites/{barcode}/{barcode}.sig_sites.annotated.gene_id.tab"
    params:
        templates_dir = get_templates_dir,
        segment = get_segment_regions,
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
        segment = get_segment_regions,
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

