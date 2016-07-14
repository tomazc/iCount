from . import annotation
from . import assemblies
from . import examples
from . import species


#wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode
# .v19.annotation.gtf.gz
#wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode
# .v19.annotation.gff3.gz


def get_genes(gtf, f_id="ID", f_parent="Parent", f_gid="gene_id",
                   f_gname="gene_name", f_gstatus="gene_status",
                   f_gtype="gene_type"):

    # "gene" records only
    genes = pybedtools.BedTool(gtf).filter(lambda r: r[2] == 'gene').saveas()

    r[2] = r.attr['gene_name']
    merged_genes = genes.merge(s=True, d=0, c=6, o='distinct').saveas()
    return genes
