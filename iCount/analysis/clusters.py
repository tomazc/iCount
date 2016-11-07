"""
Cluster sites
-------------

Merge adjacent cross-linked sites into clusters.

Read bedGraph with (significant) cross-linked sites. Cluster together sites that
are apart at most a specified number of nucleotides. Return BED file with
clusters' coordinates.
"""
import logging
import os

import pybedtools

import iCount

LOGGER = logging.getLogger(__name__)


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


def run(sites, clusters, dist=20):  # , extend=0):
    """
    Join neighboring cross-linked sites into clusters.

    Score of cluster is the sum of all its element scores.

    Parameters
    ----------
    sites : str
        Path to input BED6 file with sites.
    clusters : str
        Path to output BED6 file with merged sites.
    dist : int
        Distance between two cross_links to still merge them.

    Returns
    -------
    str
        BED file with clusters as elements.

    """
    iCount.log_inputs(LOGGER, level=logging.INFO)

    # It is required to pre-sort your data:
    sites = pybedtools.BedTool(sites).sort().saveas()

    LOGGER.info('Merging cross links form file %s', sites)
    merged = sites.merge(s=True, d=dist, c=[5, 4], o='sum,distinct').saveas()
    out = merged.sort().each(_fix_proper_bed6_format).saveas(clusters)

    LOGGER.info('Done. Results saved to: %s', os.path.abspath(out.fn))
    return out.fn
