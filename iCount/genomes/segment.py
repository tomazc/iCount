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


def get_genes(gtf_in, gtf_out, attribute='gene_name'):
    """
    Extract largest possible gene segments from input gtf file.

    Each gene can have multiple entries: for exons, introns, UTR...
    We wish to get a "maximal frame" of it's coordinates for given gene
    name, chromosome and strand.

    :param str gtf_in: absolute path to gtf input file
    :param str gtf_out: absolute path to gtf output file
    :return: sorted largest possible gene segments
    :rtype: pybedtools.BedTool
    """

    data = {}

    for interval in pybedtools.BedTool(gtf_in):
        gene_name = interval.attrs[attribute]
        chromosome = interval.chrom
        strand = interval.strand
        # Generate unique slug, since same gene name can appear on
        # multiple chromosomes or on oppsite strands.
        gene_unique_slug = '-'.join([gene_name, chromosome, strand])

        # TODO: Numerically optimze this procedure

        if gene_unique_slug in data:
            if interval.start < data[gene_unique_slug][1]:
                data[gene_unique_slug][1] = interval.start

            if interval.stop > data[gene_unique_slug][2]:
                data[gene_unique_slug][2] = interval.stop

        else:
            data[gene_unique_slug] = [interval.chrom, interval.start,
                                      interval.stop, interval.name,
                                      interval.strand]

    gs = pybedtools.BedTool(pybedtools.create_interval_from_list(
        [chrom, start, end, name, '1', strand])
        for chrom, start, end, name, strand in data.values()).saveas()

    return gs.sort().saveas(gtf_out)
