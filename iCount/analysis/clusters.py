""".. Line to protect from pydocstyle D205, D400.

Cluster sites
-------------

Merge adjacent peaks into clusters and sum cross-links within clusters.

"""
import logging
import os

import pybedtools

import iCount

LOGGER = logging.getLogger(__name__)


def _strip_empty_names(name):
    """
    Remove feature names.

    Feature name of sites are not needed when merging clusters (which should be named) and sites.
    """
    if name:
        names_list = [w.strip() for w in name.split(',')]
        names_list = [w for w in names_list if w and w != '.']
        name = ','.join(names_list)
    if name == '':
        name = '.'
    return name


def _fix_bed6_emptyname(feature):
    """
    Take a feature and convert it to BED6 format.

    http://bedtools.readthedocs.io/en/latest/content/general-usage.html
    """
    chrom = feature.chrom
    start = feature.start
    end = feature.stop
    name = feature.name
    score = feature.score
    strand = feature.strand
    name = _strip_empty_names(name)
    return pybedtools.create_interval_from_list(
        [chrom, start, end, name, score, strand])


def _fix_bed6_zeroscore_emptyname(feature):
    """
    Take a feature and convert it to BED6 format. Set score to zero.

    http://bedtools.readthedocs.io/en/latest/content/general-usage.html
    """
    chrom = feature.chrom
    start = feature.start
    end = feature.stop
    name = feature.name
    score = 0
    strand = feature.score
    name = _strip_empty_names(name)
    return pybedtools.create_interval_from_list(
        [chrom, start, end, name, score, strand])


def _select_bed6_noname(feature):
    """
    Take a feature and convert it to BED6 format. Set name to empty.

    http://bedtools.readthedocs.io/en/latest/content/general-usage.html
    """
    chrom = feature.chrom
    start = feature.start
    end = feature.stop
    name = ''
    score = feature.score
    strand = feature.strand
    name = _strip_empty_names(name)
    return pybedtools.create_interval_from_list(
        [chrom, start, end, name, score, strand])


def run(sites, peaks, clusters, dist=20, slop=3):
    """
    Join neighboring peaks (at distance dist) into clusters.

    Report sum of sites' scores within each cluster, including slop.

    Parameters
    ----------
    sites : str
        Path to input BED6 file with sites.
    peaks : str
        Path to input BED6 file with peaks (or clusters).
    clusters : str
        Path to output BED6 file with merged peaks (clusters).
    dist : int
        Distance between two peaks to merge into same cluster.
    slop : int
        Distance between site and cluster to assign site to cluster.

    Returns
    -------
    str
        BED file with clusters as elements.

    """
    iCount.log_inputs(LOGGER, level=logging.INFO)
    metrics = iCount.Metrics()

    if slop >= dist:
        LOGGER.warning('Distance between peaks (%s) should be larger than cluster slop ('
                       '%s)', dist, slop)

    # It is required to pre-sort your data:
    LOGGER.info('Reading individual sites from %s', sites)
    bt_sites = pybedtools.BedTool(sites).sort().saveas()
    LOGGER.info('Reading peaks from %s', peaks)
    bt_peaks = pybedtools.BedTool(peaks).sort().saveas()

    LOGGER.info('Merging peaks to form clusters')
    merged = bt_peaks.merge(s=True, d=dist, c=(4, 6), o='distinct,distinct').saveas()
    bt_merged = merged.each(_fix_bed6_zeroscore_emptyname).sort().saveas()

    LOGGER.info('Summing sites within identified clusters')
    # for each site, find closest cluster to which assign the site to
    bt_selected_sites = bt_sites.closest(bt_merged, s=True, d=True, t='first', stream=True).\
                                 filter(lambda b: 0 <= int(b.fields[-1]) <= slop).saveas()
    bt_selected_sites = bt_selected_sites.each(_select_bed6_noname).sort().saveas()

    # merge selected sites and previously identified clusters
    merged = bt_selected_sites.cat(bt_merged, postmerge=True,
                                   s=True, d=slop, c=[4, 5, 6], o='distinct,sum,distinct').saveas()
    out = merged.each(_fix_bed6_emptyname).sort().saveas(clusters).saveas()

    LOGGER.info('Done. Results saved to: %s', os.path.abspath(out.fn))
    return metrics
