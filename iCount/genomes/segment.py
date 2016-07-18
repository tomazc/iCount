"""
Segment genome
--------------

Parse genome annotation, segment it and prepare a number of versions needed
for mapping and in various analyses:

- regions of genes (all isoforms and other parts merged into one region)
- regions of individual region types (segment each gene into exonic,
intronic, nc, utr, etc..)
- landmarks (determine positions of exon-intron, intron-exon, exon-exon,
and other types of genomic regions)

"""

import pybedtools


# description and parameters needed for the analysis
analysis_name = 'segment'
analysis_description_short = 'segment genome regions'
analysis_description = 'Segment genome into non-overlapping regions.'

# no parameters needed
params_opt = [
    (
        'feature', 'string', 'gene', False,
        'Feature type name that specifies gene records in GTF.'
    ),
    (
        'attribute', 'string', 'gene_name', False,
        'Tag (in the attribute field of the GTF) to use for grouping records '
        'on same gene.'
    ),
]

params_pos = [
    (
        'annotation', 'GTF', 'in',
        '(input) GTF file with genome annotation.'
    ),
    (
        'segmentation_genes', 'GTF', 'out',
        '(output) GTF file with gene segments.'
    ),
]


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


def get_genes(gtf_in, gtf_out_genes, f_type='gene', a_gene_name='gene_name'):
    """Return gene coordinates in proper BED6 format.

    Extract only regions that describe genes.
    Merge overlapping regions and record associated gene names.
    All other annotation is skipped.

    """
    gs = pybedtools.BedTool(gtf_in).filter(lambda r: r[2] == f_type).saveas()
    gs1 = gs.each(_rename_to_gene_name, a_gene_name=a_gene_name).saveas()
    gs2 = gs1.sort().saveas()
    gs3 = gs2.merge(s=True, d=0, c='4', o='distinct').saveas()
    gs4 = gs3.each(_fix_proper_bed6_format).saveas(gtf_out_genes)
    return gs4
