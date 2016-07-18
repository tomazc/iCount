import pybedtools

from . import ensembl


def _rename_to_gene_name(feature, a_gene_name='gene_name'):
    chrom = feature.chrom
    start = feature.start
    end = feature.stop
    name = feature.attrs[a_gene_name]
    score = '1'
    strand = feature.strand
    # use BED6 format, see:
    # http://bedtools.readthedocs.io/en/latest/content/general-usage.html
    return pybedtools.create_interval_from_list(
        [chrom, start, end, name, score, strand]
    )


def _fix_proper_bed6_format(feature):
    chrom = feature.chrom
    start = feature.start
    end = feature.stop
    name = feature.score
    score = '1'
    strand = feature.name
    # use BED6 format, see:
    # http://bedtools.readthedocs.io/en/latest/content/general-usage.html
    return pybedtools.create_interval_from_list(
        [chrom, start, end, name, score, strand]
    )


def get_genes(gtf, f_type='gene', a_gene_name='gene_name'):
    """Return gene coordinates in proper BED6 format.

    Extract regions that describe genes.
    Merge overlapping regions and record associated gene names.
    All other annotation is skipped.

    """
    genes0 = pybedtools.BedTool(gtf).filter(lambda r: r[2] == f_type).saveas()
    genes1 = genes0.each(_rename_to_gene_name, a_gene_name=a_gene_name).saveas()
    genes2 = genes1.sort().saveas()
    genes3 = genes2.merge(s=True, d=0, c='4', o='distinct').saveas()
    genes = genes3.each(_fix_proper_bed6_format).saveas()
    return genes


def get(release, species):
    return ensembl.download_annotation(release, species)
