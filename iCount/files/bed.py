import pybedtools


def load(fname):
    return pybedtools.BedTool(fname)


def _convert_legacy_bed_format(feature):
    # old iCount legacy format:
    # chrome, start, end, [+-]score
    # where +/- indicate strand of cross-link site and score indicates the
    # intensity of interaction
    # use BED6 format, see:
    # http://bedtools.readthedocs.io/en/latest/content/general-usage.html
    chrom = feature.chrom
    start = feature.start
    end = feature.stop
    name = '.'
    if feature.name[0] == '-' or feature.name[0] == '+':
        score = feature.name[1:]
        strand = feature.name[0]
    else:
        score = feature.name
        strand = '+'
    return pybedtools.create_interval_from_list(
        [chrom, start, end, name, score, strand]
    )


def convert_legacy(bedGraph_legacy, bed_converted):
    """Convert legacy iCount's four-column format into proper BED6 format.

    Old iCount legacy format: chrome, start, end, [+-]value
    Strand can be either '+' or '-', and value indicates the intensity of
    interaction.

    The returned BED file follows the BED6 format, as explained in the
    [bedtools manual](http://bedtools.readthedocs.io/en/latest/content
    /general-usage.html).

    """
    sites = pybedtools.BedTool(bedGraph_legacy).sort().saveas()
    sites1 = sites.each(_convert_legacy_bed_format).saveas(bed_converted)
    return sites1
