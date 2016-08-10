"""
Cluster finding
---------------

Read bedGraph with (significant) cross-linked sites.
Cluster together sites that are apart at most a specified number of
nucleotides.
Return BED file with clusters' coordinates.

"""

import gzip
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
        'sites', 'BED6', 'in',
        '(input) BED6 file with cross-linked sites.'
    ),
    (
        'clusters', 'BED6', 'out',
        '(output) BED6 file with coordinates of identified clusters.'
    ),
]


def _fix_proper_bed6_format(feature):
    chrom = feature.chrom
    start = feature.start
    end = feature.stop
    name = feature.strand
    score = feature.score
    strand = feature.name
    # use BED6 format, see:
    # http://bedtools.readthedocs.io/en/latest/content/general-usage.html
    return pybedtools.create_interval_from_list(
        [chrom, start, end, name, score, strand]
    )


def run(fin_sites, fout_clusters, dist=20): #, extend=0):
    """Join neighboring cross-linked sites into clusters.

    """
    sites = pybedtools.BedTool(fin_sites).sort().saveas()
#    if extend:
#        sites = sites.slop(b=extend).saveas()
    m1 = sites.merge(s=True, d=dist, c=[5, 4], o='sum,distinct').saveas()
    m2 = m1.sort().saveas()
    m3 = m2.each(_fix_proper_bed6_format).saveas(fout_clusters)
    return m3

#     # _lowFDR.bed
#     desc = "peaks in %s, with %s permutations, %s nt neighborhood, FDR < %g, regions %s" % (

#
#     # clustered - cluster neighboring peaks, within +-15nt, report value (cDNA, tc,...) sum for region covered by cluster
#     desc = "peaks in %s, with %s permutations, FDR < %g, regions %s, clustered within %s nt neighborhood" % (

#
#     # extended - extend clusters with neighboring non-lowFDR peaks, join clusters if 15nt or less apart
#     desc = "peaks in %s, with %s permutations, FDR < %g, regions %s, clustered within %s nt neighborhood, extended to neighboring high FDR peaks or clusters %s nt or less apart" % (

