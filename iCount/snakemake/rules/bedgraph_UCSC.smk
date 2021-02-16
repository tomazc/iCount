#==============================================================================#
#     Create crosslinks bedgraphs to import on UCSC genome browser
#==============================================================================#

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
