import os

import pybedtools

from collections import Counter

from pybedtools import create_interval_from_list


# data required for CLI interface:
analysis_name = 'summary'
analysis_description_short = 'Summary report'
analysis_description = 'Report proportion of cross-link events/sites on each region type'
params_opt = [
    ('annotation_file', 'string', None, True, 'Annotation file.'),
    ('cross_links_file', 'string', None, True, 'Cross-link file (BED6).'),
    ('out_file', 'string', None, True, 'Output filename.'),
    ('types_length_file', 'string', None, False, 'File with lengths for each type.'),
]
params_pos = []


def make_types_length_file(annotation_file, out_file=None):
    """
    Calculate cumulative length for each "type" in annotation file

    In the context of this report "type" equals to the combination of 3rd
    column and attribute "biotype" from annotation file (GTF).

    :param str annotation_file: path to annotation file (should be GTF and include biotype attribute)
    :param str out_file: path to output file (if None it is determined automatically)

    :return: absolute path to out_file
    :return: str
    """
    annotation = pybedtools.BedTool(annotation_file)

    # Ignore transcript/gene lines since they are subdivided...
    excluded_types = ['transcript', 'gene']
    type_lengths = {'{} {}'.format(i[2], i.attrs.get('biotype', '.')): 0
                    for i in annotation if i[2] not in excluded_types}

    if out_file is None:
        out_file = annotation_file[:-3] + 'types_length.txt'

    for interval in annotation:
        if interval[2] not in excluded_types:
            type_ = '{} {}'.format(interval[2], interval.attrs.get('biotype', '.'))
            type_lengths[type_] += interval.length

    # Write results to file:
    with open(out_file, 'wt') as outfile:
        for type_, length in type_lengths.items():
            outfile.write('{}\t{}\n'.format(type_, length))

    return os.path.abspath(out_file)


def make_summary_report(annotation_file, cross_links_file, out_file,
                        types_length_file=None):
    """
    Make summary report from cross_links and annotation data.

    In the context of this report "type" equals to the combination of 3rd
    column and attribute "biotype" from annotation file (GTF).
    "Regions" are parts of transcript (UTR, CDS, introns...) that "non-
    intersectingly" span the whole transcript. Intergenic regions are also
    considered as region. Each region has one and only one type.

    For each of such types, number of cross-link events/sites is counted.
    Sites = sum of *all* cross-link sites that are in regions of some type.
    Events = sum of col5 for *all* cross-link sites that are in regions of
    some type (multiple events can happen on each cross link). Col5 is
    numerical value in 5th column of cross-link file.

    :param string annotation_file: path to annotation file (should be GTF and include biotype attribute)
    :param string cross_links_file: path to cross_links_file (should be BED6)
    :param string out_file: path to output file
    :returns: path to summary report file (should be equal to out_file parameter)
    :rtype: string
    """

    sites = pybedtools.BedTool(cross_links_file)
    annotation = pybedtools.BedTool(annotation_file)

    # Ignore transcript/gene lines since they are subdivided...
    excluded_types = ['transcript', 'gene']
    types_count = {'{} {}'.format(i[2], i.attrs.get('biotype', '.')): [0, 0]
                   for i in annotation if i[2] not in excluded_types}

    # Find or make file with cumulative length of types:
    if not types_length_file or not os.path.isfile(types_length_file):
        types_length_file = make_types_length_file(annotation_file)

    # read the file to dict, named type_lengths
    type_lengths = {}
    with open(types_length_file) as tfile:
        for line in tfile:
            parts = line.strip().split()
            type_lengths[' '.join(parts[:-1])] = int(parts[-1])

    overlaps = annotation.intersect(sites, sorted=True, s=True, wo=True).saveas()
    try:
        # this will raise TypeError if overlaps is empty:
        overlaps[0]
    except TypeError as error:
        raise ValueError('No intersections found. This may be caused by '
                         'different naming of chromosomes in annotation and'
                         'cross_links file ("chr1" vs. "1")')

    def finalize(temp):
        for type_, events in temp:
            types_count[type_][0] += 1 / len(temp)
            types_count[type_][1] += events / len(temp)

    temp = []
    previous_site = None

    for overlap in overlaps:
        if overlap[2] not in excluded_types:
            site = overlap[10]  # take start of cross-link as unique indentifier for site
            if site != previous_site:
                finalize(temp)
                temp = []
            type_ = '{} {}'.format(overlap[2], overlap.attrs.get('biotype', '.'))
            temp.append([type_, int(overlap[13])])
            site = previous_site

    finalize(temp)

    # Produce report file:
    header = ['type', 'length', 'length %', 'sites #', 'sites %',
              'sites enrichment', 'events #', 'events %', 'events enrichment']
    sum_sites = sum([i[0] for i in types_count.values()])
    sum_events = sum([i[1] for i in types_count.values()])
    total_length = sum([i for i in type_lengths.values()])
    with open(out_file, 'wt') as out:
        out.write('\t'.join(header) + '\n')
        for type_, [sites, events] in sorted(types_count.items()):
            length_percent = type_lengths[type_] / total_length
            # site_percent = '{0:.3f}'.format(sites / sum_sites)
            site_percent = sites / sum_sites
            site_enrichment = site_percent / length_percent
            event_percent = events / sum_events
            event_enrichment = event_percent / length_percent
            line = [type_, type_lengths[type_], length_percent,
                    sites, site_percent, site_enrichment,
                    events, event_percent, event_enrichment]
            out.write('\t'.join(map(str, line)) + '\n')

    return os.path.abspath(out_file)
