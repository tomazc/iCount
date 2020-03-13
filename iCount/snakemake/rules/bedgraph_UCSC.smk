#==============================================================================#
#     Create crosslinks bedgraphs to import on UCSC genome browser
#==============================================================================#

def bedgraph_header(wildcards):
    # db=\"{mapto}\" removed
    # return ("{project}_{sample_name}_{protein}_{method}_{mapto}".format(project=config['project'], sample_name=samples.loc[wildcards.barcode, "sample_name"], mapto=samples.loc[wildcards.barcode, "mapto"], method=samples.loc[wildcards.barcode, "method"],	protein=samples.loc[wildcards.barcode, "protein"],	cells_tissue=samples.loc[wildcards.barcode, "cells_tissue"],	condition=samples.loc[wildcards.barcode, "condition"],))
    return ("track type=bedGraph name=\"{project}_{sample_name}_{protein}_{method}_{mapto}_unique.xl.bedgraph.bed\" description=\"{project}_{sample_name}_{protein}_{method}_{mapto}\" "
            "color=\"120,101,172\"  mapped_to=\"{mapto}\" altColor=\"200,120,59\" lib_id=\"{project}\" maxHeightPixels=\"100:50:0\" visibility=\"full\" "
            "tissue=\"{cells_tissue}\" protein=\"{protein}\" species=\"{mapto}\" condition=\"{condition}\" res_type=\"T\" priority=\"20\" \n".format(project=config['project'], sample_name=samples.loc[wildcards.barcode, "sample_name"], mapto=samples.loc[wildcards.barcode, "mapto"], method=samples.loc[wildcards.barcode, "method"],	protein=samples.loc[wildcards.barcode, "protein"],	cells_tissue=samples.loc[wildcards.barcode, "cells_tissue"],	condition=samples.loc[wildcards.barcode, "condition"],))


rule bedgraphUCSC:
    input:
        xlsites="{project}/xlsites/{barcode}/{barcode}.unique.xl.bed",
    output:
        bedgraph="{project}/xlsites/{barcode}/{barcode}.unique.xl.UCSC.bedgraph",
    params:
        header=bedgraph_header,
        genome=get_genome,
    run:
        # Convert ENSMBL to UCSC. Thanks to Devon Ryan for creation of mapping tables
        # If your genome is not included please check: https://github.com/dpryan79/ChromosomeMappings
        d = {}
        if params.genome == 'homo_sapiens':
            f = open("data/GRCh38_ensembl2UCSC.txt")
        elif params.genome == 'mus_musculus':
            f = open("data/GRCm38_ensembl2UCSC.txt")
        else:
            print("Please add mapping file to trasnform your genome coordinates to UCSC compatible chromosomes")
            f = ""

        for line in f:
            cols = line.strip().split("\t")
            if len(cols) < 2 or cols[1] == "":
                continue
            d[cols[0]] = cols[1]

        f.close()

        fin = open(input.xlsites)
        fout = open(output.bedgraph, "w")
        fout.write(params.header)
        line = fin.readline()
        while line:
            col = line.rstrip('\n').rsplit('\t')
            chr = col[0]
            count = col[4]
            strand = col[5]
            if strand == '-':
                count = '-' + count
            fout.write(d[chr] + '\t' + col[1] + '\t' + col[2] + '\t' + count + '\n')
            line = fin.readline()

        fin.close()
        fout.close()
