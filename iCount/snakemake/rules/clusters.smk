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


