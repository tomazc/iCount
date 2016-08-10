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


import pybedtools

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


def _get_non_cds_exons(cdses, exons, intervals):
    """
    Identify (parts of) exons that are not CDS and mark them as UTR

    UTRs upstream/downstream start_codon/stop_codon are named
    "UTR5"/"UTR3", respectively.

    Two scenarions have to be indentified for each exon:

        * no CDS is overlapping exon - exon is outside start-stop codon region
        * part of exon overlaps with CDS, but part does not

    :param list cdses: list of all CDS in transcript group
    :param list exons: list of all exons in transcript group
    :param list intervals: list of all intervals in transcript group
    :return: list of UTR intervals in transcript group
    :rtype: list
    """
    utrs = []
    int0 = intervals[0]
    strand = int0.strand
    start_codon = next((i for i in intervals if i.fields[2] == 'start_codon'), None)
    stop_codon = next((i for i in intervals if i.fields[2] == 'stop_codon'), None)
    if not start_codon:
        if strand == '+':
            start = min([c.start for c in cdses])
            stop = start + 3
        else:
            stop = max([c.stop for c in cdses])
            start = stop - 3
        start_codon = create_interval_from_list(
            int0[:2] + ['start_codon', start, stop, '.', strand, '.', int0[8]])

    if not stop_codon:
        if strand == '+':
            start = max([c.stop for c in cdses])
            stop = start + 3
        else:
            stop = min([c.start for c in cdses])
            start = stop - 3
        stop_codon = create_interval_from_list(
            int0[:2] + ['stop_codon', start, stop, '.', strand, '.', int0[8]])

    for exon in exons:
        # Identify the cds that is inside given exon:
        cds = next((c for c in cdses if _a_in_b(c, exon)), None)
        if not cds:
            # This exon is UTR completely -just determine if UTR3 or UTR5
            if (strand == '+' and exon.start < start_codon.start) or \
               (strand == '-' and exon.stop > start_codon.stop):
                mode = 'UTR5'
            else:
                mode = 'UTR3'
            start, stop = exon.start + 1, exon.stop
            utrs.append(create_interval_from_list(
                exon[:2] + [mode, start, stop, '.', strand, '.', exon[8]]))

        elif not (abs(cds.start - exon.start) <= 3 and abs(cds.stop - exon.stop) <= 3):
            # CDS and exon partially overlap - determine borders and UTR5/UTR3:
            if strand == '+' and exon.stop > stop_codon.stop:
                    mode = 'UTR3'
                    start = stop_codon.stop + 1
                    stop = exon.stop
            elif strand == '+' and exon.start < start_codon.start:
                    mode = 'UTR5'
                    start = exon.start + 1
                    stop = start_codon.start
            elif strand == '-' and exon.start < stop_codon.start:
                    mode = 'UTR3'
                    start = exon.start + 1
                    stop = stop_codon.start
            elif strand == '-' and exon.stop > start_codon.stop:
                    mode = 'UTR5'
                    start = start_codon.stop + 1
                    stop = exon.stop
            else:
                print(cds)
                print(exon)
                print(strand)
                for i in intervals:
                    print(i)
                raise ValueError('No strand info or something strange...')

            utrs.append(create_interval_from_list(
                exon[:2] + [mode, start, stop, '.', strand, '.', exon[8]]))
        else:
            # exon and CDS overlap completely or just for stop codon apart
            # This case is already covered in _process_transcript_group
            pass

    return utrs


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
    start_stop = [(e1.stop + 1, e2.start) for e1, e2 in zip(exons, exons[1:])]

    # For introns, keep only a subset of key-value pairs from column 8:
    intron_keys = ['gene_id', 'gene_name', 'transcript_id', 'transcript_name']
    col8 = ' '.join(['{} "{}";'.format(key, value) for key, value in
                    exons[0].attrs.items() if key in intron_keys])

    ex1 = exons[0]
    strand = ex1.strand
    return [create_interval_from_list(
            ex1[:2] + ['intron', i[0], i[1], '.', strand, '.', col8])
            for i in start_stop]


def _process_transcript_group(intervals):
    """
    Get a list of interval objects representing the transcript group

    These objects respect the "convention", defined in get_regions docstring.

    :param list intervals:
    :return: Modified list of pybedtools interval objects
    :rtype: list
    """
    # Container for interval objects
    regions = []

    # Extract the interval describing the whole transcript
    index = next(i for i in range(len(intervals)) if intervals[i][2] == 'transcript')
    regions.append(intervals.pop(index))

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
        cdses = [create_interval_from_list(
            i[:2] + ['CDS'] + i[3:]) for i in intervals if i[2] == 'CDS']
        # check that all CDSs are within exons:
        for cds in cdses:
            assert(any([_a_in_b(cds, exon) for exon in exons]))

        # Include CDS regions:
        regions.extend(cdses)
        # Include UTR regions:
        regions.extend(_get_non_cds_exons(cdses, exons, intervals))

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
    col8 = 'gene_id "."; transcript_id ".";'
    gtf = pybedtools.BedTool(create_interval_from_list(
        [i[0], '.', 'intergenic', str(int(i[1]) + 1), i[2], '.', strand, '.', col8])
        for i in intergenic_bed).saveas()

    return os.path.abspath(gtf.fn)


def get_regions(gtf_in, gtf_out, genome_file):
    """
    Create new gtf file with custom annotation and filtered content.

    Each line in new file should define one of the following elements:

        * gene
        * transcript
        * CDS
        * intron
        * UTR3
        * UTR5
        * ncRNA
        * intergenic

    Name of third field (interval.fields[2]) should correspond to one of theese names.

    :param string gtf_in: path to input GTF file
    :param string gtf_in: path to output GTF file
    :return: absolute path to output GTF file
    :rtype: str
    """
    # Containers for grouping of transcripts and for output data
    transcript_group = []
    data = []
    # List to keep track of all already processed transcript id's:
    transcript_ids = []
    no_transcript_id = 0

    # Only consider the chromosomes that are present in genome file
    chromosomes = [line.strip().split()[0] for line in open(genome_file)]

    length = pybedtools.BedTool(gtf_in).count()
    j = 0

    # for interval in [i for i in pybedtools.BedTool(gtf_in) if i.chrom in chromosomes]:
    for interval in pybedtools.BedTool(gtf_in):
        j += 1
        sys.stdout.write("\r{0:.1f} %".format(j / length * 100))
        sys.stdout.flush()
        # Only check for the chromosomes provided in genome file:
        if interval.chrom in chromosomes:
            # Segments without 'transcript_id' attributes are the ones that
            # define genes - checked for homo_sapiens.84 ENSEMBL release
            if 'transcript_id' not in interval.attrs:
                no_transcript_id += 1
                data.append(interval)
                # NOTE: If desired, only lines that define genes can be obtained
                # from output file by: grep -v "transcript_id" gtf_out.gtf
            else:
                if not transcript_group or interval.attrs['transcript_id'] == transcript_ids[-1]:
                    # Interval belongs to the current transcript group
                    transcript_group.append(interval)
                    if not transcript_ids:
                        transcript_ids.append(interval.attrs['transcript_id'])
                else:
                    # This is the first interval of new transcript!
                    # But first process all the data for previous transcript:
                    data.extend(_process_transcript_group(transcript_group))

                    # Same transcript should be written in one single block.
                    # This ensures that when new transcript_id is found, one
                    # can be sure to have all the items from previous transcript
                    # group. Report if this is not the case:
                    if interval.attrs['transcript_id'] in transcript_ids:
                        raise ValueError("Transcripts are not ordered!")
                    transcript_ids.append(interval.attrs['transcript_id'])

                    # Make list for new group of transcript items and add first one:
                    transcript_group = [interval]

    print('There are {} intervals without "transcript_id" attribute.'.format(
        no_transcript_id))

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
