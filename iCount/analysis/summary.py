"""
Cross-link site summary
-----------------------

Report proportion of cross-link events/sites on each region type
"""
import os
import re
import logging

import pybedtools

import iCount

from pybedtools import create_interval_from_list

LOGGER = logging.getLogger(__name__)


def make_types_length_file(annotation_file, out_file=None, subtype='biotype',
                           excluded_types=None):
    """
    Calculate the number of non-overlapping base pairs of each "type".

    In the context of this function "type" equals to the combination of 3rd
    column and attribute subtype from annotation file (GTF).


    Parameters
    ----------
    annotation_file : str
        Path to annotation file (should be GTF and include subtype attribute).
    out_file : str
        Path to output file (if None it is determined automatically).
    subtype : int
        Subtype.
    excluded_types : list_str
        Excluded types.

    Returns
    -------
    str
        Absolute path to out_file

    """
    excluded_types = excluded_types or []
    annotation = pybedtools.BedTool(annotation_file).filter(
        lambda x: x[2] not in excluded_types).sort().saveas()

    if out_file is None:
        match = re.match(r'([\w_]+\.\d+.).*', os.path.basename(annotation_file))
        if match:
            out_file = str(match.group(1)) + 'types_length.txt'
        else:
            out_file = annotation_file + 'types_length.txt'

    data = {}
    for interval in annotation:
        if subtype:
            type_ = '{} {}'.format(interval[2], interval.attrs.get(subtype, '.'))
        else:
            type_ = interval[2]
        if type_ not in data:
            data[type_] = []
        data[type_].append(create_interval_from_list(interval.fields))

    type_lengths = {}
    for type_, list_ in data.items():
        total_bp = 0
        for feature in pybedtools.BedTool(i for i in list_).merge(s=True):
            total_bp += len(feature)
        type_lengths[type_] = total_bp

    # Write results to file:
    with open(out_file, 'wt') as outfile:
        for type_, length in sorted(type_lengths.items()):
            outfile.write('{}\t{}\n'.format(type_, length))

    return os.path.abspath(out_file)


def make_summary_report(annotation_file, cross_links_file, out_file,
                        chrom_length_file, types_length_file=None, ndigits='8',
                        subtype='biotype',
                        excluded_types=None):
    """
    Make summary report from cross-link and annotation data.

    In the context of this report "type" equals to the combination of 3rd
    column and attribute `subtype` from annotation file (GTF).
    "Regions" are parts of transcript (UTR, CDS, introns...) that "non-
    intersectingly" span the whole transcript. Intergenic regions are also
    considered as region. Each region has one and only one type.

    For each of such types, number of cross-link events/sites is counted.
    Sites = sum of *all* cross-link sites that are in regions of some type.
    Events = sum of col5 for *all* cross-link sites that are in regions of
    some type (multiple events can happen on each cross link). Col5 is
    numerical value in 5th column of cross-link file.

    Parameters
    ----------
    annotation_file : str
        Path to annotation file (should be GTF and include subtype attribute).
    cross_links_file : str
        Path to cross_links_file (should be BED6)
    out_file : str
        Path to output file
    chrom_length_file : str
        Path to file with chromosome lengths
    types_length_file : str
        Path to file with lengths of each type
    ndigits : int
        Number of decimal places in results
    subtype : str
        Name of attribute to be used as subtype
    excluded_types : list_str
        Types in 3rd column to exclude form analysis

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

    # If not given/present, make file with cumulative length for each type:
    if not types_length_file or not os.path.isfile(types_length_file):
        LOGGER.info('types_length_file not given - calculating it')
        types_length_file = make_types_length_file(
            annotation_file, subtype=subtype, excluded_types=excluded_types)

    # read the file to dict, named type_lengths
    type_lengths = {}
    with open(types_length_file) as tfile:
        for line in tfile:
            parts = line.strip().split()
            type_lengths[' '.join(parts[:-1])] = int(parts[-1])

    # sorted=True - invokes memory efficient algorithm for large files
    # s=True - only report hits in B that overlap A on the same strand
    # wb=True - Write the original entry in B for each overlap
    LOGGER.info('Calculating intersection between cross-link and annotation_file...')
    overlaps = cross_links.intersect(annotation, sorted=True, s=True, wb=True).saveas()
    try:
        # this will raise TypeError if overlaps is empty:
        overlaps[0]
    except (IndexError, TypeError):
        raise ValueError('No intersections found. This may be caused by '
                         'different naming of chromosomes in annotation and '
                         'cross-links file (example: "chr1" vs. "1")')

    # dict structure = type_: [# of sites, # of events]
    type_counter = {type_: [0, 0] for type_ in type_lengths.keys()}
    site_types = []
    previous_segment = overlaps[0]

    def finalize(types, segment):
        """Increase counter for all types that intersect with segment"""
        for type_ in set(types):
            type_counter[type_][0] += 1  # sites
            type_counter[type_][1] += int(segment[4])  # events

    LOGGER.info('Extracting summary from data...')
    for segment in overlaps:
        # detect if segment contains new cross-link site:
        if segment.start != previous_segment.start or segment.strand != previous_segment.strand:
            finalize(site_types, previous_segment)
            site_types = []
        if subtype:
            # Extract subtype attribute:
            stype = re.match(r'.*{} "(.*)";'.format(subtype), segment[-1])
            site_types.append('{} {}'.format(segment[8], stype.group(1) if stype else '.'))
        else:
            site_types.append(segment[8])
        previous_segment = segment
    finalize(site_types, previous_segment)

    # Produce report file:
    header = ['type', 'length', 'length %', 'sites #', 'sites %',
              'sites enrichment', 'events #', 'events %', 'events enrichment']
    sum_sites = sum([i[0] for i in type_counter.values()])
    sum_events = sum([i[1] for i in type_counter.values()])

    # total genome len = sum_of_all_chrom_length * 2 (there are + and - strand):
    total_length = sum([int(line.strip().split()[1]) for line in
                        open(chrom_length_file)]) * 2
    with open(out_file, 'wt') as out:
        out.write('\t'.join(header) + '\n')
        for type_, [sites, events] in sorted(type_counter.items()):
            length_percent = type_lengths[type_] / total_length
            site_percent = sites / sum_sites
            site_enrichment = site_percent / length_percent
            event_percent = events / sum_events
            event_enrichment = event_percent / length_percent
            line = [type_, type_lengths[type_], length_percent,
                    sites, site_percent, site_enrichment,
                    events, event_percent, event_enrichment]
            line = line[:1] + [round(i, int(ndigits)) for i in line[1:]]
            out.write('\t'.join(map(str, line)) + '\n')

    LOGGER.info('Done. Results written in %s.', os.path.abspath(out_file))
    return os.path.abspath(out_file)
