"""
Cluster finding
---------------

Read bedGraph with (significant) cross-linked sites.
Cluster together sites that are apart at most a specified number of
nucleotides.
Return BED file with clusters' coordinates.

"""

import pybedtools

# description and parameters needed for the analysis
analysis_name = 'clusters'
analysis_description_short = 'cluster analysis'
analysis_description = 'Merge adjacent cross-linked sites into clusters.'

params_opt = [
    (
        'dist', 'int_range', (20, 0, 200), False,
        'Maximum distance between sites that can be merged into a single '
        'cluster.'
    ),
# TODO: pybedtools requires info on chromosome sizes
#    (   'extend', 'int_range', (0, 0, 200), False,
#        'Before merging, extend both upstream and downstream of actual sites '
#        'by a specified number of nucleotides.'
#    ),
]

params_pos = [
    (
        'sites', 'BED6', 'in', 1,
        '(input) BED6 file with cross-linked sites.'
    ),
    (
        'clusters', 'BED6', 'out', 1,
        '(output) BED6 file with coordinates of identified clusters.'
    ),
]


def _fix_proper_bed6_format(feature):
    """
    Take a feature and convert it to BED6 format

    http://bedtools.readthedocs.io/en/latest/content/general-usage.html
    """
    chrom = feature.chrom
    start = feature.start
    end = feature.stop
    name = feature.strand
    score = feature.score
    strand = feature.name
    return pybedtools.create_interval_from_list(
        [chrom, start, end, name, score, strand])


def run(fin_sites, fout_clusters, dist=20):  # , extend=0):
    """
    Join neighboring cross-linked sites into clusters.

    Score of cluster is the sum of all its element scores.

    :param str fin_sites: Path to input file with sites (BED6 format)
    :param str fout_clusters: Path to output file with merged sites (BED6)
    :param int dist: Distance between two cross_links to still merge them.
    :return: BED file with clusters as elements
    :rtype: str
    """
    # It is required to pre-sort your data:
    sites = pybedtools.BedTool(fin_sites).sort().saveas()
    # if extend:
    #    sites = sites.slop(b=extend).saveas()
    merged = sites.merge(s=True, d=dist, c=[5, 4], o='sum,distinct').saveas()
    out = merged.sort().each(_fix_proper_bed6_format).saveas(fout_clusters)
    return out.fn

#     # _lowFDR.bed
#     desc = "peaks in %s, with %s permutations, %s nt neighborhood, FDR < %g, regions %s" % (

#
#     # clustered - cluster neighboring peaks, within +-15nt, report value (cDNA, tc,...) sum for region covered by cluster
#     desc = "peaks in %s, with %s permutations, FDR < %g, regions %s, clustered within %s nt neighborhood" % (

#
#     # extended - extend clusters with neighboring non-lowFDR peaks, join clusters if 15nt or less apart
#     desc = "peaks in %s, with %s permutations, FDR < %g, regions %s, clustered within %s nt neighborhood, extended to neighboring high FDR peaks or clusters %s nt or less apart" % (

