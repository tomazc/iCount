"""
Cross-link site annotation
--------------------------

Annotate each cross link site with types of regions that intersect with it.

"""
import re
import os
import logging

import pybedtools

import iCount

from pybedtools import create_interval_from_list

LOGGER = logging.getLogger(__name__)


def annotate_cross_links(annotation_file, cross_links_file, out_file,
                         subtype='biotype',
                         excluded_types=None):
    """
    Annotate each cross-link site with all region types that intersect it.

    Each cross link can overlap/intersect with many intervals in annotation
    file. Each of these intervals has a given type. Make a copy of cross
    link file and include all these types in 4th column.

    In the context of this report "type" equals to the combination of 3rd column
    and attribute `subtype` from annotation file (GTF). Regions/intervals are
    parts of transcript (UTR, CDS, introns...) that "non- intersectingly" span
    the whole transcript. However, since transcripts can overlap, also intervals
    belonging to different transcripts can overlap. Intergenic regions are also
    considered as region. Each region has one and only one type.

    Parameters
    ----------
    annotation_file : str
        Path to annotation file (should be GTF and include `subtype` attribute).
    cross_links_file : str
        Path to cross_links_file (should be BED6).
    out_file : str
        Path to output file.
    subtype : str
        Subtype.
    excluded_types : list_str
        Excluded types.

    Returns
    -------
    str
        Path to summary report file (should be equal to out_file parameter)

    """
    iCount.log_inputs(LOGGER, level=logging.INFO)

    excluded_types = excluded_types or []
    cross_links = pybedtools.BedTool(cross_links_file).sort().saveas()
    annotation = pybedtools.BedTool(annotation_file).filter(
        lambda x: x[2] not in excluded_types).sort().saveas()

    LOGGER.info('Calculating overlaps between cross-link and annotation_file...')
    overlaps = cross_links.intersect(annotation, sorted=True, s=True, wb=True).saveas()
    try:
        # this will raise TypeError if overlaps is empty:
        overlaps[0]
    except (IndexError, TypeError):
        raise ValueError('No intersections found. This may be caused by '
                         'different naming of chromosomes in annotation and'
                         'cross_links file ("chr1" vs. "1")')

    data = []  # cotainer for final annotated BED file intervals
    site_types = []  # cotainer for all types intersecting with given cross-link
    previous_interval = overlaps[0]

    def finalize(types, site):
        """Make annotated (with all intersecting types) cross link interval."""
        data.append(create_interval_from_list(
            site[0:3] + [', '.join(map(str, sorted(set(types))))] + site[4:6]))

    for interval in overlaps:
        # Detect new cross link:
        if interval.start != previous_interval.start or interval.strand != previous_interval.strand:
            finalize(site_types, previous_interval)
            site_types = []
        if subtype:
            stype = re.match(r'.*{} "(.*?)";'.format(subtype), interval[-1])  # Extract subtype attribute
            site_types.append('{} {}'.format(interval[8], stype.group(1) if stype else '.'))
        else:
            site_types.append(interval[8])
        previous_interval = interval
    finalize(site_types, previous_interval)

    # Produce annotated cross-link file:
    LOGGER.info('Writing results to file...')
    pybedtools.BedTool(line for line in data).saveas(out_file)
    LOGGER.info('Done. Output saved to: %s', os.path.abspath(out_file))
    return os.path.abspath(out_file)
