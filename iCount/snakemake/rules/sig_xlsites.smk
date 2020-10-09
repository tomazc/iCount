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

