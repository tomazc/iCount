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

import os
import sys
import time
import shutil
import tempfile
import multiprocessing

import pybedtools

from collections import Counter

from pybedtools import create_interval_from_list

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


def _a_in_b(a, b):
    """
    Check if interval a is inside interval b

    :param pybedtools.Interval a: interval
    :param pybedtools.Interval b: interval
    :return: true if a is inside b, false otherwise
    :rtype: bool
    """
    return a.start >= b.start and a.stop <= b.stop


def _add_biotype_attribute(gene_content_):
    """
    Add `biotype` attribute to all intervals in gene_content

    biotype attribute is equal to transcript_biotype value if present,
    else gene_biotype if present else value in column 2 (index 1). The
    last option can happen in some of early releases.

    :param dict gene_content_: intervals in gene separated by transcript id
    :return: same gene_content_ object with added `biotype` attributes
    :rtype: dict
    """
    gene_content = gene_content_.copy()

    # Determine gene biotype:
    gbiotype = gene_content['gene'].attrs.get('gene_biotype', None)

    # List to keep track of all possible biotypes in gene:
    gene_biotypes = [gbiotype] if gbiotype else []

    for transcript_id, transcript_intervals in gene_content.items():
        if transcript_id == 'gene':
            continue
        new_intervals = []

        exon = [i for i in transcript_intervals if i[2] in ['CDS', 'ncRNA']][0]
        gbiotype = exon.attrs.get('gene_biotype', None)
        tbiotype = exon.attrs.get('transcript_biotype', None)
        biotype = tbiotype if tbiotype else (gbiotype if gbiotype else exon[1])
        gene_biotypes.append(biotype)
        for interval in transcript_intervals:
            new_intervals.append(create_interval_from_list(
                interval[:8] + [interval[8] + ' biotype "{}";'.format(biotype)]))
        gene_content[transcript_id] = new_intervals

    # Finally, make also gene biotype: a list of all biotypes in
    # gene, sorted by frequency:
    biotype = ', '.join([i[0] for i in sorted(
        Counter(gene_biotypes).items(), key=lambda x: x[1], reverse=True)])
    interval = gene_content['gene']
    gene_content['gene'] = create_interval_from_list(
        interval[:8] + [interval[8] + ' biotype "[{}]";'.format(biotype)])

    return gene_content


def _check_consistency(intervals):
    """
    Check that intervals in transcript are consistent with theory.

    The test checks two things:

        * interval[i].stop == interval[i + 1].start - this ensures that
          intervals are covering whole transcript and do not overlap
        * that type of given interval is consistent with the type of
          the preceding one

    :param list intervals: intervals representing transcript
    :return: None
    :rtype: None
    :raises: AssertionError: if any of the assert statements fails
    """
    can_follow = {
        '+': {
            'UTR5': ['intron', 'CDS'],
            'CDS': ['intron', 'UTR3'],
            'intron': ['CDS', 'ncRNA', 'UTR3', 'UTR5'],
            'UTR3': ['intron'],
            'ncRNA': ['intron']},
        '-': {
            'UTR3': ['intron', 'CDS'],
            'CDS': ['intron', 'UTR5'],
            'intron': ['CDS', 'ncRNA', 'UTR3', 'UTR5'],
            'UTR5': ['intron'],
            'ncRNA': ['intron']}}
    intervals = intervals.copy()
    index = next(i for i in range(len(intervals)) if intervals[i][2] == 'transcript')
    transcript_interval = intervals.pop(index)
    strand = intervals[0].strand

    intervals = sorted(intervals, key=lambda x: x.start)
    assert(transcript_interval.start == intervals[0].start)
    assert(transcript_interval.stop == intervals[-1].stop)
    for first, second in zip(intervals, intervals[1:]):
        assert(first.stop == second.start)
        assert(second[2] in can_follow[strand][first[2]])


def _get_non_cds_exons(cdses, exons, intervals):
    """
    Identify (parts of) exons that are not CDS and classify them

    The intervals in the output list of this function should have
    "UTR3" or "UTR5" in their third field.

    First, all stop codons and CDS regions are merged where posssible.
    For each stop codon there two options:
        * touching CDS - merge them
        * NOT touching CDS - rename to CDS
        * stop_codon and CDS overlap completely - ignore stop codon

    For *each* exon, many cases are posssible. Here is a "tree" of all
    possible cases for positive strand:

    * CDS inside exon
        * CDS and exon overlap perfectly  #1 - no UTR here
        * CDS and exon do NOT overlap perfectly
            * cds.start != exon.start & cds.stop == exon.stop
                * exon = utr5 + cds  #2
            * cds.start == exon.start & cds.stop != exon.stop
                * exon = cds + utr3  #3
            * cds.start != exon.start & cds.stop != exon.stop
                * exon = utr5 + cds + utr3  #4
    * CDS not inside exon
        * exon completely UTR, just determine wheather UTR3 or UTR5  #5

    Currently, all 5 cases are covered.

    :param list cdses: list of all CDS in transcript group
    :param list exons: list of all exons in transcript group
    :param list intervals: list of all intervals in transcript group
    :return: list of UTR intervals in transcript group
    :rtype: list
    """
    utrs = []
    int0 = intervals[0]
    strand = int0.strand
    stop_codons = [i for i in intervals if i.fields[2] == 'stop_codon']
    cdses = cdses.copy()

    # Merge stop_codons with cds where posssible:
    replace_cds_indexes = []
    replace_cdses = []
    new_cdses = []
    for stop_codon in stop_codons:
        touching_cds = next((
            cds for cds in enumerate(cdses) if
            stop_codon.stop == cds[1].start or
            cds[1].stop == stop_codon.start), None)
        # in some GTF files stop_codon is also cds, identify such cases:
        complete_overlap = any(cds.start == stop_codon.start and cds.stop == stop_codon.stop for cds in cdses)
        if complete_overlap:
            continue  # in this case stop_codon is already turned in cds...
        elif touching_cds:
            cds_index, cds = touching_cds
            replace_cds_indexes.append(cds_index)
            start = min(cds.start, stop_codon.start) + 1
            stop = max(cds.stop, stop_codon.stop)
            replace_cdses.append(create_interval_from_list(
                cds[:3] + [start, stop] + cds[5:]))
        else:
            new_cdses.append(create_interval_from_list(
                stop_codon[:2] + ['CDS'] + stop_codon[3:]))
    for index, cds in zip(replace_cds_indexes, replace_cdses):
        cdses[index] = cds
    cdses.extend(new_cdses)

    for exon in exons:
        if not any(_a_in_b(cds, exon) for cds in cdses):
            # no CDS in exon - completely UTR! Just determine UTR3/UTR5:
            if ((exon.strand == '+' and exon.stop >= max([i.stop for i in cdses])) or
                    (exon.strand == '-' and exon.start <= min([i.start for i in cdses]))):
                mode = "UTR3"
            else:
                mode = "UTR5"
            utrs.append(create_interval_from_list(
                exon[:2] + [mode, exon.start + 1, exon.stop, '.', strand, '.', exon[8]]))

        else:
            # CDS in exons! Identify which one:
            cds = next((c for c in cdses if _a_in_b(c, exon)), None)
            assert(cds is not None)

            if cds.start != exon.start:
                # UTR in the beggining:
                mode = 'UTR5' if exon.strand == '+' else "UTR3"
                utrs.append(create_interval_from_list(
                    exon[:2] + [mode, exon.start + 1, cds.start, '.', strand, '.', exon[8]]))
            if cds.stop != exon.stop:
                # UTR in the end:
                mode = "UTR3" if exon.strand == '+' else "UTR5"
                utrs.append(create_interval_from_list(
                    exon[:2] + [mode, cds.stop + 1, exon.stop, '.', strand, '.', exon[8]]))

    return cdses, utrs


def filter_col8(interval, keys=None):
    """
    Filter the content of last column in interval

    :param pybedtools.Interval interval: interval
    :param list keys: list of keys that should be kept in result
    :return: filtered content of last column
    :rtype: string
    """
    if keys is None:
        keys = ['gene_id', 'gene_name', 'transcript_id', 'transcript_name']
    return ' '.join(['{} "{}";'.format(key, value) for key, value in
                    interval.attrs.items() if key in keys])


def _get_introns(exons):
    """
    Calculate the positions of introns and return them as a list of intervals.

    Introns are regions located between exons. When constructing them,
    one can copy all the data from exons except:

        * start and stop codon have to be determined
        * score and frame column are set to undefined value: '.'
        * the name of third column is set to 'intron'
        * content referring to exon is removed from last column

    :param list exons: list of exons, sorted by coordinates
    :return: list of pybedtools interval objects representing introns
    :rtype: list
    """
    # start of intron is on exon1.stop + 1
    # stop of intron is on exon2.start (since in GTF: e2.start = e2[3] - 1)
    start_stop = [(e1.stop + 1, e2.start) for e1, e2 in zip(exons, exons[1:])]

    # For introns, keep only a subset of key-value pairs from column 8:
    col8 = filter_col8(exons[0])

    ex1 = exons[0]
    strand = ex1.strand
    return [create_interval_from_list(
            ex1[:2] + ['intron', i[0], i[1], '.', strand, '.', col8])
            for i in start_stop]


def _process_transcript_group(intervals):
    """
    Process a list of intervals in the transcript group

    These intervals should not overlap and should completely span the
    transcript. Their third column should be one of the following:

        * transcript
        * CDS
        * intron
        * UTR3
        * UTR5
        * ncRNA

    :param list intervals:
    :return: Modified list of pybedtools interval objects
    :rtype: list
    """
    # Container for interval objects
    regions = []

    # Get the interval describing whole transcript (make it if not present)
    index = next((i for i in range(len(intervals)) if intervals[i][2] == 'transcript'), None)
    if index is not None:
        regions.append(intervals.pop(index))
    else:
        # Manually create "transcript interval":
        i1 = intervals[0]
        col8 = filter_col8(i1)

        start = min([i.start for i in intervals])
        stop = max([i.stop for i in intervals])
        regions.append(create_interval_from_list(
            i1[:2] + ['transcript', start + 1, stop] + i1[5:8] + [col8]))

    exons = [i for i in intervals if i[2] == 'exon']
    assert(len(exons) != 0)

    # Sort exones by exom number (reverse if strand == '-'):
    exons = sorted(exons, key=lambda x: int(x.attrs['exon_number']),
                   reverse=exons[0].strand == '-')
    # Confirm that they are really sorted: get their starts:
    starts = [exon.start for exon in exons]
    assert(starts == sorted(starts))

    # Gaps between exons are introns. Apend introns to regions container:
    regions.extend(_get_introns(exons))

    if not {'CDS', 'start_codon', 'stop_codon'} & {i[2] for i in intervals}:
        # If no CDS/stop_codon/start_codon transcript name should be ncRNA.
        regions.extend([create_interval_from_list(i[:2] + ['ncRNA'] + i[3:]) for i in intervals])
    else:
        cdses = [i for i in intervals if i[2] == 'CDS']
        # check that all CDSs are within exons:
        for cds in cdses:
            assert(any([_a_in_b(cds, exon) for exon in exons]))

        # Determine UTR regions and new_cds (cds joined with stop
        # codons where possible):
        new_cdses, utrs = _get_non_cds_exons(cdses, exons, intervals)
        regions.extend(utrs)
        regions.extend(new_cdses)

    # check the consistency of regions and make a report if check fails:
    try:
        _check_consistency(regions)
    except AssertionError as e:
        print('X' * 80)
        for i in sorted(intervals, key=lambda x: x.start):
            print(i)
        print('*' * 80)
        for i in sorted(regions, key=lambda x: x.start):
            print(i)
        raise e

    return regions


def _complement(gs, genome_file, strand):
    """
    Get the complement of intervals in gs that have strand == `strand`

    Required structure of genome_file: first column has to be chromosome
    name and the second column has to be chromosome length. Files produced
    with `samtools faidx` (*.fai file extension) respect this rule.

    Possible options for strand param: '+', '-' and '.'.

    :param str gs: path to GTF file with content
    :param str genome_file: path to genome_file (*.fai or similar)
    :param string strand: Strand for which to compute complement.
    :return: absolute path to GTF file with complement segments.
    :rtype: str
    """
    assert(strand in ['+', '-', '.'])

    # Filter the content file by strand
    gs_strand_only = pybedtools.BedTool(gs).filter(
        lambda r: r.strand == strand).saveas()

    # Sorted input is required for complement calculation:
    gs_sorted = gs_strand_only.sort(faidx=genome_file).saveas()

    # This step will fail if gs_sorted has different set of chromosomes
    # as they are defined in genome_file:
    intergenic_bed = gs_sorted.complement(g=genome_file).saveas()

    # intergenic_bed is BED3 file. We need to make it GTF:
    # Note the differences in BED and  GTF:
    # https://pythonhosted.org/pybedtools/intervals.html#bed-is-0-based-others-are-1-based
    # This effectively means:
    # gtf_start = bed_start + 1
    # gtf_stop = bed_stop
    col8 = 'gene_id "."; transcript_id ".";'
    gtf = pybedtools.BedTool(create_interval_from_list(
        [i[0], '.', 'intergenic', str(int(i[1]) + 1), i[2], '.', strand, '.', col8])
        for i in intergenic_bed).saveas()

    return os.path.abspath(gtf.fn)


def _get_gene_content(gtf_in, chromosomes, show_progress=False):
    """
    Generator giving groups of intervals belonging to one gene

    The yielded structure in each iteration is a dictionary that has
    key-value pairs:

        * 'gene': interval if type gene
        * 'transcript_id#1': intervals corresponding to transcript_id#1
        * 'transcript_id#2': intervals corresponding to transcript_id#2
        ...

    :param str gtf_in: path to gtf input file
    :param list chromosomes: list of chromosomes to consider
    :param bool show_progress: switch to show progress
    :return: all intervals in gene, separated by transcript_id
    :rtype: dict
    """
    # Lists to keep track of all already processed genes/transcripts:
    gene_ids = []
    transcript_ids = []

    current_transcript = None
    current_gene = None
    gene_content = {}

    length = pybedtools.BedTool(gtf_in).count()
    j = 0

    def finalize(gene_content):
        """
        Procedure before returning group of intervals belonging to one gene
        """
        if 'gene' not in gene_content:
            # Manually create "gene interval":
            i1 = next(iter(gene_content.values()))[0]
            col8 = filter_col8(i1)
            start = min([i.start for j in gene_content.values() for i in j])
            stop = max([i.stop for j in gene_content.values() for i in j])
            gene_content['gene'] = create_interval_from_list(
                i1[:2] + ['gene', start + 1, stop] + i1[5:8] + [col8])
        return gene_content

    previous_progress = ''
    for interval in pybedtools.BedTool(gtf_in):
        j += 1
        if show_progress:
            progress = '\r{0:.1f} %'.format(j / length * 100)
            if progress != previous_progress:
                print(progress, end="", flush=True)
                previous_progress = progress

        if interval.chrom in chromosomes:
            # Segments without 'transcript_id' attributes are the ones that
            # define genes. such intervals are not in all releases.
            if interval.attrs['gene_id'] == current_gene:
                if interval.attrs['transcript_id'] == current_transcript:
                    # Same gene, same transcript: just add to container:
                    gene_content[current_transcript].append(interval)
                else:
                    # New transcript - confirm that it is really a new one:
                    current_transcript = interval.attrs['transcript_id']
                    assert(current_transcript not in transcript_ids)
                    transcript_ids.append(current_transcript)
                    gene_content[current_transcript] = [interval]

            else:  # New gene!
                # First process old content:
                if gene_content:  # To survive the first iteration
                    yield finalize(gene_content)

                # Confirm that it is really new gene!
                assert(interval.attrs['gene_id'] not in gene_ids)
                # Then add it to already processed genes:
                current_gene = interval.attrs['gene_id']
                gene_ids.append(current_gene)

                # Make empty container and classify interval
                gene_content = {}
                if interval[2] == 'gene':
                    gene_content['gene'] = interval
                elif 'transcript_id' in interval.attrs:
                    current_transcript = interval.attrs['transcript_id']
                    assert(current_transcript not in transcript_ids)
                    transcript_ids.append(current_transcript)
                    gene_content[current_transcript] = [interval]
                else:
                    raise Exception("Unexpected situation!")

    # for the last iteration:
    yield finalize(gene_content)


def get_regions(gtf_in, gtf_out, genome_file, cores=1, show_progress=False):
    """
    Create new gtf file with custom annotation and filtered content.

    Each line in new file should define one of the following elements:

        * gene
        * transcript
        * CDS
        * intron
        * UTR3
        * UTR5
        * stop_codon
        * ncRNA
        * intergenic

    Name of third field (interval.fields[2]) should correspond to one
    of theese names. Only consider GTF entries of chromosomes given in
    genome_file.

    :param string gtf_in: path to input GTF file
    :param string gtf_out: path to output GTF file
    :param string genome_file: path to genome_file (*.fai file or similar )
    :param int cores: number of computer cpu-s to use for calculation
    :param bool show_progress: switch to show progress
    :return: absolute path to output GTF file
    :rtype: str
    """
    # Container for storing intermediate data
    data = []

    with open(genome_file) as gfile:
        chromosomes = [line.strip().split()[0] for line in gfile]

    def process_gene(gene_content):
        """
        Process each group of intervals belonging to gene

        Process each transcript_group in gene_content, add 'biotype'
        attribute to all intervals and include them in `data`.
        """
        assert('gene' in gene_content)

        for id_, transcript_group in gene_content.items():
            if id_ == 'gene':
                continue
            gene_content[id_] = _process_transcript_group(transcript_group)

        # Add biotype attribute to all intervals:
        gene_content = _add_biotype_attribute(gene_content)

        for id_, transcript_group in gene_content.items():
            if id_ == 'gene':
                continue
            data.extend(transcript_group)
        data.append(gene_content['gene'])

    for gene_content in _get_gene_content(gtf_in, chromosomes, show_progress):
        process_gene(gene_content)
    # This can be replaced with: multiprocessing.Pool, but it causes huge
    # memory usage. Possible explanation and solution:
    # http://stackoverflow.com/questions/21485319/high-memory-usage-using-python-multiprocessing
    # p = multiprocessing.Pool(cores, maxtasksperchild=100)
    # p.map(process_gene, _get_gene_content(gtf_in, chromosomes))

    # Produce GTF/GFF file from data:
    gs = pybedtools.BedTool(i.fields for i in data).saveas()

    intergenic_pos = _complement(gs.fn, genome_file, '+')
    intergenic_neg = _complement(gs.fn, genome_file, '-')

    # Join the gs, intergenic_pos and intergenic_neg in one file:
    f2 = tempfile.NamedTemporaryFile(delete=False)
    for infile in [gs.fn, intergenic_pos, intergenic_neg]:
            shutil.copyfileobj(open(infile, 'rb'), f2)
    f2.close()

    return pybedtools.BedTool(f2.name).sort().saveas(gtf_out).fn
