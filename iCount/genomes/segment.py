""".. Line to protect from pydocstyle D205, D400.

Segmentation
------------

Parse annotation file into internal iCount structure - segmentation.

Currently, only annotations from ENSEMBl and GENCODE are supported.
http://www.gencodegenes.org/
http://www.ensembl.org

Segmentation is used in almost all further analyses.

In segmentation, each transcript is partitioned into so called
regions/intervals. Such regions must span the whole transcript, but should not
intersect with each other. However, higher hierarchy levels: transcripts and
genes can of course intersect each other.

Example of possible segmentation::

    Genome level: |---------------------------------------------------|

    Gene level:    |--------------gene1--------------|   |-intergenic-|
                                 |---------gene2--------|

    Transcript l.: |----------transcript1---------|
                           |-------transcript2-------|
                                 |------transcript3-----|

    Region level:  |-CDS-||-intron-||-CDS-||-UTR3-|

For simplicity, only the partition of transcript1 is presented.

"""
import os
import shutil
import logging
import tempfile

from collections import Counter

import pybedtools

from pybedtools import create_interval_from_list

import iCount

LOGGER = logging.getLogger(__name__)


def _first_two_columns(input_file):
    """Keep just first two columns of file."""
    ofile = tempfile.NamedTemporaryFile(mode='wt', delete=False)
    with open(input_file) as ifile:
        for line in ifile:
            col_1_2 = line.strip().split()[:2]
            ofile.write('\t'.join(col_1_2) + '\n')
    ofile.close()
    return os.path.abspath(ofile.name)


def _a_in_b(first, second):
    """
    Check if interval a is inside interval b.

    Parameters
    ----------
    a : pybedtools.Interval
        First interval.
    b : pybedtools.Interval
        Second interval.

    Returns
    -------
    bool
        True if a is inside b, false otherwise.
    """
    return first.start >= second.start and first.stop <= second.stop


def _get_biotype(interval):
    """
    Get interval biotype.

    Biotype of interval is equal to transcript biotype value if present,
    else gene biotype if present else value in column 2 (index 1). The
    last option can only happen in some of early ENSEMBL releases.

    Transcript_biotype and gene biotype values are determined from attributes:

        * In ENSEMBL annotations, they are stored in
          ``transcript_biotype`` and ``transcript_type`` attributes.
        * In GENCODE annotations, they are stored in
          ``gene_biotype`` and ``gene_type`` attributes.

    Parameters
    ----------
    interval : pybedtools.Interval
        Interval == line in GTf file.

    Returns
    -------
    string
        Biotype of interval.

    """
    if "transcript_biotype" in interval.attrs:
        return interval.attrs['transcript_biotype']
    elif "transcript_type" in interval.attrs:
        return interval.attrs['transcript_type']
    elif "gene_biotype" in interval.attrs:
        return interval.attrs['gene_biotype']
    elif "gene_type" in interval.attrs:
        return interval.attrs['gene_type']
    else:
        return interval[1]


def _add_biotype_value(interval, biotype):
    """Add biotype value to interval."""
    col8 = interval[8] if interval[8] != '.' else ''
    return create_interval_from_list(
        interval[:8] + [col8 + ' biotype "{}";'.format(biotype)])


def _add_biotype_attribute(gene_content):
    """
    Add `biotype` attribute to all intervals in gene_content.

    Parameters
    ----------
    gene_content_ : dict
        Intervals in gene separated by transcript id.

    Returns
    -------
    dict
        Same gene_content_ object with added `biotype` attributes.
    """
    gene_content = gene_content.copy()

    # Determine gene biotype:
    gbiotype = _get_biotype(gene_content['gene'])
    # List to keep track of all possible biotypes in gene:
    gene_biotypes = [gbiotype] if gbiotype else []

    for transcript_id, transcript_intervals in gene_content.items():
        if transcript_id == 'gene':
            continue
        first_exon = [i for i in transcript_intervals if i[2] in ['CDS', 'ncRNA']][0]
        biotype = _get_biotype(first_exon)
        gene_biotypes.append(biotype)

        new_intervals = []
        for interval in transcript_intervals:
            new_intervals.append(_add_biotype_value(interval, biotype))
        gene_content[transcript_id] = new_intervals

    # Finally, make also gene biotype: a list of all biotypes in gene,
    # sorted by frequency. Additionally, another sorting is added to sort
    # by alphabet if counts are equal.
    biotype = ', '.join([i[0] for i in sorted(
        sorted(Counter(gene_biotypes).items()), key=lambda x: x[1], reverse=True)])
    gene_content['gene'] = _add_biotype_value(gene_content['gene'], biotype)

    return gene_content


def _check_consistency(intervals):
    """
    Check that intervals in transcript are consistent with theory.

    The test checks two things:

        * interval[i].stop == interval[i + 1].start - this ensures that
          intervals are covering whole transcript and do not overlap
        * that type of given interval is consistent with the type of
          the preceding one


    Parameters
    ----------
    intervals : list
        Intervals representing transcript.

    Returns
    -------
    None
        None.

    Raises
    ------
    AssertionError
        IN case any of the assert statements fails.

    """
    can_follow = {
        '+': {
            'UTR5': ['intron', 'CDS'],
            'CDS': ['intron', 'UTR3'],
            'intron': ['CDS', 'ncRNA', 'UTR3', 'UTR5'],
            'UTR3': ['intron'],
            'ncRNA': ['intron'],
        },
        '-': {
            'UTR3': ['intron', 'CDS'],
            'CDS': ['intron', 'UTR5'],
            'intron': ['CDS', 'ncRNA', 'UTR3', 'UTR5'],
            'UTR5': ['intron'],
            'ncRNA': ['intron'],
        }
    }
    intervals = intervals.copy()
    try:
        index = next(i for i in range(len(intervals)) if
                     intervals[i][2] == 'transcript')
    except StopIteration:
        raise ValueError("No transcript interval in list of intervals.")
    transcript_interval = intervals.pop(index)
    strand = intervals[0].strand

    intervals = sorted(intervals, key=lambda x: x.start)
    assert transcript_interval.start == intervals[0].start
    assert transcript_interval.stop == intervals[-1].stop
    for first, second in zip(intervals, intervals[1:]):
        assert first.stop == second.start
        assert second[2] in can_follow[strand][first[2]]


def _get_non_cds_exons(cdses, exons, intervals):
    """
    Identify (parts of) exons that are not CDS and classify them.

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


    Parameters
    ----------
    cdses : list
        List of all CDS in transcript group.
    exons : list
        List of all exons in transcript group.
    intervals : list
        List of all intervals in transcript group.

    Returns
    -------
    list
        List of UTR intervals in transcript group.

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
        complete_overlap = any(
            cds.start == stop_codon.start and cds.stop == stop_codon.stop for cds in cdses)
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
            assert cds is not None

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


def _filter_col8(interval, keys=None):
    """
    Filter the content of last column in interval.

    Parameters
    ----------
    interval : pybedtools.Interval
        Interval.
    keys : list
        List of keys that should be kept in result.

    Returns
    -------
    string
        Filtered content of last column.

    """
    if keys is None:
        keys = ['gene_id', 'gene_name', 'transcript_id', 'transcript_name']
    return ' '.join(['{} "{}";'.format(key, value) for key, value in
                     sorted(interval.attrs.items()) if key in keys])


def _get_introns(exons):
    """
    Calculate the positions of introns and return them as a list of intervals.

    Introns are regions located between exons. When constructing them,
    one can copy all the data from exons except:

        * start and stop bp for each intron have to be determined
        * score and frame column are set to undefined value: '.'
        * the name of third column is set to 'intron'
        * content referring to exon is removed from last column


    Parameters
    ----------
    exons : list
        List of exons, sorted by coordinates.

    Returns
    -------
    list
        List of pybedtools interval objects representing introns.

    """
    # start of intron is on exon1.stop + 1
    # stop of intron is on exon2.start (since in GTF: e2.start = e2[3] - 1)
    start_stop = [(e1.stop + 1, e2.start) for e1, e2 in zip(exons, exons[1:])]

    # For introns, keep only a subset of key-value pairs from column 8:
    col8 = _filter_col8(exons[0])

    ex1 = exons[0]
    strand = ex1.strand
    return [create_interval_from_list(
        ex1[:2] + ['intron', i[0], i[1], '.', strand, '.', col8]) for i in start_stop]


def _process_transcript_group(intervals):
    """
    Process a list of intervals in the transcript group.

    Processed intervals should not overlap and should completely span the
    transcript. Their third column should be one of the following:

        * transcript
        * CDS
        * intron
        * UTR3
        * UTR5
        * ncRNA

    Parameters
    ----------
    intervals : list
        List of intervals in transcript to process.

    Returns
    -------
    list
        Modified list of pybedtools interval objects.
    """
    # Container for interval objects
    regions = []

    # Get the interval describing whole transcript (make it if not present)
    index = next((i for i in range(len(intervals)) if intervals[i][2] == 'transcript'), None)
    if index is not None:
        regions.append(intervals.pop(index))
    else:
        # Manually create "transcript interval":
        int1 = intervals[0]
        col8 = _filter_col8(int1)

        start = min([i.start for i in intervals])
        stop = max([i.stop for i in intervals])
        regions.append(create_interval_from_list(
            int1[:2] + ['transcript', start + 1, stop] + int1[5:8] + [col8]))

    exons = [i for i in intervals if i[2] == 'exon']
    assert len(exons) != 0

    # Sort exones by exom number (reverse if strand == '-'):
    exons = sorted(exons, key=lambda x: int(x.attrs['exon_number']),
                   reverse=exons[0].strand == '-')
    # Confirm that they are really sorted: get their starts:
    starts = [exon.start for exon in exons]
    assert starts == sorted(starts)

    # Gaps between exons are introns. Apend introns to regions container:
    regions.extend(_get_introns(exons))

    if not {'CDS', 'start_codon', 'stop_codon'} & {i[2] for i in intervals}:
        # If no CDS/stop_codon/start_codon transcript name should be ncRNA.
        regions.extend([create_interval_from_list(i[:2] + ['ncRNA'] + i[3:]) for i in intervals])
    else:
        cdses = [i for i in intervals if i[2] == 'CDS']
        # check that all CDSs are within exons:
        for cds in cdses:
            assert any([_a_in_b(cds, exon) for exon in exons])

        # Determine UTR regions and new_cds (cds joined with stop
        # codons where possible):
        new_cdses, utrs = _get_non_cds_exons(cdses, exons, intervals)
        regions.extend(utrs)
        regions.extend(new_cdses)

    # check the consistency of regions and make a report if check fails:
    try:
        _check_consistency(regions)
    except AssertionError as error:
        print('X' * 80)
        for i in sorted(intervals, key=lambda x: x.start):
            print(i)
        print('*' * 80)
        for i in sorted(regions, key=lambda x: x.start):
            print(i)
        raise error

    return regions


def _complement(gtf, genome_file, strand, type_name='intergenic'):
    """
    Get the complement of intervals in gtf that have strand == `strand`.

    Required structure of genome_file: first column has to be chromosome
    name and the second column has to be chromosome length. Files produced
    with `samtools faidx` (*.fai file extension) respect this rule.

    Possible options for strand param: '+', '-' and '.'.


    Parameters
    ----------
    gtf : str
        Path to GTF file with content.
    genome_file : str
        Path to genome_file (*.fai or similar).
    strand : string
        Strand for which to compute complement.

    Returns
    -------
    str
        Absolute path to GTF file with complement segments.

    """
    assert(strand in ['+', '-', '.'])

    # Filter the content file by strand
    gtf_strand_only = pybedtools.BedTool(gtf).filter(
        lambda r: r[6] == strand).saveas()

    # pylint: disable=unexpected-keyword-arg
    # Sorted input is required for complement calculation:
    gtf_sorted = gtf_strand_only.sort(faidx=genome_file).saveas()

    # This step will fail if gtf_sorted has different set of chromosomes
    # as they are defined in genome_file:
    intergenic_bed = gtf_sorted.complement(g=genome_file).saveas()

    # intergenic_bed is BED3 file. We need to make it GTF:
    # Note the differences in BED and  GTF:
    # https://pythonhosted.org/pybedtools/intervals.html#bed-is-0-based-others-are-1-based
    # This effectively means:
    # gtf_start = bed_start + 1
    # gtf_stop = bed_stop
    col8 = 'gene_id "."; transcript_id ".";'
    gtf = pybedtools.BedTool(
        create_interval_from_list(
            [i[0], '.', type_name, str(int(i[1]) + 1), i[2], '.', strand, '.', col8]
        )
        for i in intergenic_bed
    ).saveas()

    return os.path.abspath(gtf.fn)


def _get_gene_content(gtf, chromosomes, report_progress=False):
    """
    Generator giving groups of intervals belonging to one gene.

    The yielded structure in each iteration is a dictionary that has
    key-value pairs:

        * 'gene': interval if type gene
        * 'transcript_id#1': intervals corresponding to transcript_id#1
        * 'transcript_id#2': intervals corresponding to transcript_id#2
        ...


    Parameters
    ----------
    gtf : str
        Path to gtf input file.
    chromosomes : list
        List of chromosomes to consider.
    report_progress : bool
        Switch to show progress.

    Returns
    -------
    dict
        All intervals in gene, separated by transcript_id.

    """
    # Lists to keep track of all already processed genes/transcripts:
    gene_ids = []
    transcript_ids = []

    current_transcript = None
    current_gene = None
    gene_content = {}

    def finalize(gene_content):
        """Procedure before returning group of intervals belonging to one gene."""
        if 'gene' not in gene_content:
            # Manually create "gene interval":
            int1 = next(iter(gene_content.values()))[0]
            col8 = _filter_col8(int1)
            start = min([i.start for j in gene_content.values() for i in j])
            stop = max([i.stop for j in gene_content.values() for i in j])
            gene_content['gene'] = create_interval_from_list(
                int1[:2] + ['gene', start + 1, stop] + int1[5:8] + [col8])
        return gene_content

    length = pybedtools.BedTool(gtf).count()
    progress, j = 0, 0
    for interval in pybedtools.BedTool(gtf):
        j += 1
        if report_progress:
            new_progress = j / length
            # pylint: disable=protected-access
            progress = iCount._log_progress(new_progress, progress, LOGGER)

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
                    assert current_transcript not in transcript_ids
                    transcript_ids.append(current_transcript)
                    gene_content[current_transcript] = [interval]

            else:  # New gene!
                # First process old content:
                if gene_content:  # To survive the first iteration
                    yield finalize(gene_content)

                # Confirm that it is really new gene!
                assert interval.attrs['gene_id'] not in gene_ids
                # Then add it to already processed genes:
                current_gene = interval.attrs['gene_id']
                gene_ids.append(current_gene)

                # Make empty container and classify interval
                gene_content = {}
                if interval[2] == 'gene':
                    gene_content['gene'] = interval
                elif 'transcript_id' in interval.attrs:
                    current_transcript = interval.attrs['transcript_id']
                    assert current_transcript not in transcript_ids
                    transcript_ids.append(current_transcript)
                    gene_content[current_transcript] = [interval]
                else:
                    raise Exception("Unexpected situation!")

    # for the last iteration:
    yield finalize(gene_content)


def get_regions(annotation, segmentation, fai, report_progress=False):
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

    Parameters
    ----------
    annotation : str
        Path to input GTF file.
    segmentation : str
        Path to output GTF file.
    fai : str
        Path to input genome_file (.fai or similar).
    report_progress : bool
        Switch to show progress.

    Returns
    -------
    str
        Absolute path to output GTF file.

    """
    iCount.log_inputs(LOGGER, level=logging.INFO)

    # Container for storing intermediate data
    data = []

    LOGGER.debug('Opening genome file: %s', fai)

    # Keep just first two column or _complement will fail...
    fai = _first_two_columns(fai)
    with open(fai) as gfile:
        chromosomes = [line.strip().split()[0] for line in gfile]

    def process_gene(gene_content):
        """
        Process each group of intervals belonging to gene.

        Process each transcript_group in gene_content, add 'biotype'
        attribute to all intervals and include them in `data`.
        """
        assert 'gene' in gene_content

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

    LOGGER.debug('Processing genome annotation from: %s', annotation)
    for gene_content in _get_gene_content(annotation, chromosomes, report_progress):
        process_gene(gene_content)
        LOGGER.debug('Just processed gene: %s', gene_content['gene'].attrs['gene_id'])
    # This can be replaced with: multiprocessing.Pool, but it causes huge
    # memory usage. Possible explanation and solution:
    # http://stackoverflow.com/questions/21485319/high-memory-usage-using-python-multiprocessing
    # TODO: check and fix execution in parallel, look example
    # https://daler.github.io/pybedtools/3-brief-examples.html#example-3-count-reads-in-introns-and-exons-in-parallel
    # p = multiprocessing.Pool(threads, maxtasksperchild=100)
    # p.map(process_gene, _get_gene_content(gtf_in, chromosomes))

    # Produce GTF/GFF file from data:
    gtf = pybedtools.BedTool(i.fields for i in data).saveas()

    LOGGER.info('Calculating intergenic regions...')
    intergenic_pos = _complement(gtf.fn, fai, '+')
    intergenic_neg = _complement(gtf.fn, fai, '-')

    # Join the gtf, intergenic_pos and intergenic_neg in one file:
    file2 = tempfile.NamedTemporaryFile(delete=False)
    for infile in [gtf.fn, intergenic_pos, intergenic_neg]:
        shutil.copyfileobj(open(infile, 'rb'), file2)
    file2.close()

    file3 = pybedtools.BedTool(file2.name).sort().saveas(segmentation)
    LOGGER.info('Segmentation stored in %s', file3.fn)
    return file3.fn


def _prepare_annotation(ann_file):
    """
    Parse annotation file to hierarchical structure.

    Utility function to transformn internal annotation file to the following
    hierarchical structure::

        annotation = {
            (chr1, +): {
                gene_id#1: {
                    'gene_segment': gene_segment,
                    transcript_id#1: [transcript_segment, exon1, intron1, exon2, ...],
                    transcript_id#2: [transcript_segment, exon1, intron1, exon2, ...],
                    ...
                },
                gene_id#2: {},
                ...
            },
            (chr1, -),
            (chr2, +),
            ...
        }

    Note that intergenic segments have multiple roles: the same intergenic
    segment has the role of gene, transcript and sub-transcript segment. This
    eases the treatment of intergenic regions in algorithms that use this
    function (rnamaps, xlsites, ...)

    Parameters
    ----------
    ann_file : str
        Path to GTF file, produces by ``get_regions`` function.

    Returns
    -------
    dict
        Annotation, wrapped in dict with chrom-strand/gene/transcript levels of
        depth.
    """
    annotation = {}

    for segment in pybedtools.BedTool(ann_file):
        if segment[2] == 'gene':
            annotation.setdefault((segment.chrom, segment.strand), {}). \
                setdefault(segment.attrs['gene_id'], {}). \
                setdefault('gene_segment', segment)

        elif segment[2] == 'intergenic':
            # Make artificial_id from chromosome, strand and start:
            fake_gid = 'G_{}_{}_{}'.format(segment.chrom, segment.strand, segment.start)
            fake_tid = 'T_{}_{}_{}'.format(segment.chrom, segment.strand, segment.start)
            annotation.setdefault((segment.chrom, segment.strand), {}). \
                setdefault(fake_gid, {})['gene_segment'] = segment
            annotation[(segment.chrom, segment.strand)][fake_gid][fake_tid] = [segment]

        elif segment[2] == 'transcript':
            annotation.setdefault((segment.chrom, segment.strand), {}). \
                setdefault(segment.attrs['gene_id'], {}). \
                setdefault(segment.attrs['transcript_id'], []). \
                insert(0, segment)  # Ensure that transcript segment is the first one in list.

        else:  # normal segment
            annotation.setdefault((segment.chrom, segment.strand), {}). \
                setdefault(segment.attrs['gene_id'], {}). \
                setdefault(segment.attrs['transcript_id'], []). \
                append(segment)

    return annotation
