import os
import re

import pybedtools

from pybedtools import create_interval_from_list


# data required for CLI interface:
analysis_name = 'summary'
analysis_description_short = 'Summary report'
analysis_description = 'Report proportion of cross-link events/sites on each region type'
params_opt = [
    ('annotation_file', 'string', None, True, 'Annotation file.'),
    ('cross_links_file', 'string', None, True, 'Cross-link file (BED6).'),
    ('out_file', 'string', None, True, 'Output filename.'),
    ('chrom_length_file', 'string', None, True, 'File with chromosome lengths.'),
    ('types_length_file', 'string', None, False, 'File with lengths for each type.'),
    ('ndigits', 'string', '8', False, 'Number of decimal places in results'),
    ('subtype', 'string', 'biotype', False, 'Attribute defining subtype'),
    ('excluded_types', 'str_list', [], False, 'Types from third column to exclude from analysis'),
]
params_pos = []


def make_types_length_file(annotation_file, out_file=None, subtype='biotype',
                           excluded_types=None):
    """
    Calculate the number of non-overlapping base pairs of each "type".

    In the context of this function "type" equals to the combination of 3rd
    column and attribute subtype from annotation file (GTF).

    :param str annotation_file: path to annotation file (should be GTF and
    include subtype attribute)
    :param str out_file: path to output file (if None it is determined automatically)
    :return: absolute path to out_file
    :return: str
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

    :param string annotation_file: path to annotation file (should be GTF and
    include subtype attribute)
    :param string cross_links_file: path to cross_links_file (should be BED6)
    :param string out_file: path to output file
    :param string chrom_length_file: path to file with chromosome lengths
    :param string types_length_file: path to file with lengths of each type
    :param string ndigits: Number of decimal places in results
    :param string subtype: name of attribute to be used as subtype
    :param list excluded_types: types in 3rd column to exclude form analysis
    :returns: path to summary report file (should be equal to out_file parameter)
    :rtype: string
    """
    excluded_types = excluded_types or []
    cross_links = pybedtools.BedTool(cross_links_file).sort().saveas()
    annotation = pybedtools.BedTool(annotation_file).filter(
        lambda x: x[2] not in excluded_types).sort().saveas()

    # If not given/present, make file with cumulative length for each type:
    if not types_length_file or not os.path.isfile(types_length_file):
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

    for segment in overlaps:
        # detect if segment contains new cross-link site:
        if segment.start != previous_segment.start or segment.strand != previous_segment.strand:
            finalize(site_types, previous_segment)
            site_types = []
        if subtype:
            stype = re.match(r'.*{} "(.*)";'.format(subtype), segment[-1])  # Extract subtype attribute
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

    return os.path.abspath(out_file)
