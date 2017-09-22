""".. Line to protect from pydocstyle D205, D400.

Identify and quantify cross-linked sites
----------------------------------------

Quantity cross-link events and determine their positions.

The simplest version of this script would oprate on such example::

    |--a---b--- reference sequence, chr 14, positive strand ------------
        |rbc1---R1--------|
        |rbc1---R2------|
        |rbc2---R3------|
            |rbc1--------R4-------------|
            |rbc3------R5-----|

Five reads (R1-R5) are mapped to a reference sequence (chromosome 14, positive
strand). Reads start on two distinct positons. On first position, there is
R1-R3. Cross-link site is located one nucleotide before start of the read (on
negative strand, one nucleotide after end of read). However, we wish to count
number of cDNA molecules, not the number of reads. This can be done by counting
the number of distinct random barcodes (sometimes also called randomers). So in
upper example, we have:

    Postion a: 3 reads, 2 distinct random barcodes = 2 cDNA's
    Postion b: 2 reads, 2 distinct random barcodes = 2 cDNA's

However, things can get complicated when a single read is mapped in multiple
parts. This can happen for several reasons. One common example is that introns
are removed during transcription. This can be illustrated with the following
image::

    |---------------- reference -----------------------

    |--------------------------transcript--------------------------|
    |---UTR5---||---intron---||---exon---||---intron---||---exon---|

                                |-------R1-------|
                                |--R2.1-->              <-R2.2-|
                      |-R3.1->      <-R3.2->            <-R3.3-|
                        |-R4.1->      <-R4.2-|

Read R1 and R2 are starting on same position. For the sake of argument, let's
also pretend they also have same random barcode. In so, we would count them as a
single cDNA molecule (= single cross-link event), even though it is obvious that
they represent two separate cross-link events. In order to fix this, we count
not just the number of different randomers on same position, but also number of
different "second-start" coordinates. Second-start coordinate is just the
coordinate of the second part of the read. This way, the actual number of
cross-link events can be determined more accurately. If read is not split, it's
second-start coordinate is 0. If read has multiple "holes" (as read R3) we
determine second-start from the largest hole.

Reads whose second-start do NOT fall on segmentation (like R4) are stored in a
separate BAM file ``sites_strange``. They should be treated with special care,
since they can indicate not-yet annotated features in genome. If segmentation is
not given, all reads with holes, bigger than ``holesize_th`` are considered
strange.

Another parameter needs more explanation: ``group_by``. When algorithm starts,
reads from BAM file are grouped in hierarchical structure by::

    * chromosome and strand
    * cross-link position
    * random barcode
    * second-start

Each second-start group receives 1 cDNA score. This score is divided to each
read in group (if there are 5 reads in group, each one gets 1/5 score). This
enables that each read has it's cDNA score and of course, 1 "read score". This
scores can be assigned to start (actually, to cross-link position), midlle or
end position of read. By default, score is of course assigned to cross-link
location. But for diagnostic purpuses, scores can also be assigned to middle or
end coordinate of the read.



TODO: check overlap between unique and multimap BED files, should be small,
otherwise, we should think of a more approapriate quantification of (division
of randomers among) unique mapped sites that overlap with multimapped reads

"""
import re
import os
import math
import logging

import pybedtools
import pysam
from pysam import AlignmentFile  # pylint: disable=no-name-in-module

import iCount
from iCount.files import _f2s, get_temp_file_name


LOGGER = logging.getLogger(__name__)
VALID_NUCLEOTIDES = set('ATCGN')
RANDOM_BARCODE_REGEX = r'.*:rbc:([ATCGN]+).*'


def _iter_bed_dict(bed, val_index=None):
    """Iterate through dict object."""
    if val_index is not None:
        for (chrome, strand), by_pos in bed.items():
            for pos, val in by_pos.items():
                val = val[val_index]
                yield pybedtools.create_interval_from_list(
                    [chrome, pos, pos + 1, '.', _f2s(val), strand]
                )
    else:
        for (chrome, strand), by_pos in bed.items():
            for pos, val in by_pos.items():
                yield pybedtools.create_interval_from_list(
                    [chrome, pos, pos + 1, '.', _f2s(val), strand]
                )


def _save_dict(bed, out_fname, val_index=None):
    """Save data from dict to BED file."""
    sites = pybedtools.BedTool(
        _iter_bed_dict(bed, val_index=val_index)
    ).saveas()
    sites1 = sites.sort().saveas(out_fname)
    return sites1


def _get_random_barcode(query_name, metrics):
    """Extract random barcode from ``query_name``."""
    match = re.match(RANDOM_BARCODE_REGEX, query_name)
    if match:
        barcode = match.group(1)
    elif ':' in query_name:
        barcode = query_name.rsplit(':', 1)[1]
        if set(barcode) - VALID_NUCLEOTIDES:
            # invalid barcode characters
            barcode = ''
            metrics.invalidrandomer_recs += 1
    else:
        barcode = ''
        metrics.norandomer_recs += 1

    return barcode


def _match(seq1, seq2, mismatches):
    """
    Test if sequence seq1 and seq2 are sufficiently similar.

    Parameters
    ----------
    seq1 : str
        First sequence.
    seq2 : str
        Second sequence.
    mismatches : int
        Number of allowed mismatches between given sequences.

    Returns
    -------
    bool
        Do sequence `seq1` and `seq2` have less or equal than ``mismatches``

    """
    seq1, seq2 = seq1.upper(), seq2.upper()
    matches = sum([(nuc1 == 'N' or nuc2 == 'N' or nuc1 == nuc2) for nuc1, nuc2 in zip(seq1, seq2)])
    return max(len(seq1), len(seq2)) - matches <= mismatches


def _update(cur_vals, to_add):
    """
    Add the values from ``to_add`` to appropriate place in ``cur_vals``.

    Note: cur_vals is updated in place!

    Parameters
    ----------
    to_add : dict
        Dict with update data.

    Returns
    -------
    None
        None.

    """
    for pos, vals_to_add in to_add.items():
        prev_vals = cur_vals.get(pos, [0] * len(vals_to_add))
        cur_vals[pos] = [p + n for p, n in zip(prev_vals, vals_to_add)]


def _merge_similar_randomers(by_bc, mismatches, ratio_th=0.1):
    """
    Merge randomers on same site that are max ``mismatches`` different.

    Input parameter `by_bc` has te following structure:
    by_bc = {
        'AAA': [(middle_pos, end_pos, read_len, num_mapped, second_start),  # hit1
                (middle_pos, end_pos, read_len, num_mapped, second_start),  # hit2
                (middle_pos, end_pos, read_len, num_mapped, second_start),  # ...
        ]
        'AAT': [(middle_pos, end_pos, read_len, num_mapped, second_start),  # hit1
                (middle_pos, end_pos, read_len, num_mapped, second_start),  # hit2
        ]

    Steps in function:
        0. Identify ambiguous randomers ('N' characters in barcode)

        1. For each ambiguous randomer, identify similar non-ambiguous one. If
        match is found, move hits form ambiguous to the non-ambiguous one. If no
        match is found, declare ambiguous randomer (one that has 'N's) as
        non-ambiguous anyway.

        2. Identify and accept randomers with strong support.
        Randomers that are supported by a threshold number of hits are accepted as
        unique. Threshold is defined as the proportion (ratio_th) of the number of
        reads assigned to the most frequent randomer.

        3. Merge remaining.
        For each barcode with number of reads below the threshold, identify if there
        exists any similar one. If there is, join the hits form second barcode to
        the first one.

    TODO: Code should be improved in step #1. Instead of finding any match,
    match with least difference should be found. Check also the skipped unit
    test in tests/test_xlsites.py

    Parameters
    ----------
    by_bc : dict
        Dictionary of barcodes and their hits.
    mismatches : int
        Reads on same position with random barcode differing less than
        ``mismatches`` are grouped together.

    Returns
    -------
    None
        None, since input `by_bc` is modified in-place.

    """
    # Step #1: assign ambigious randomers to non-ambiguous randomers
    # Identify ambiguous barcodes first.
    nonambig_bcs = set()  # accepted_barcodes
    ambig_bcs = []  # ambiguous_barcodes
    for barcode in by_bc.keys():
        undefined_nucleotides = barcode.count('N')
        if undefined_nucleotides == 0:
            nonambig_bcs.add(barcode)
        else:
            ambig_bcs.append((undefined_nucleotides, barcode))

    # For each ambiguous randomer, identify similar non-ambiguous one.
    # If match is found, move hits form ambiguous to the non-ambiguous one.
    # If no match is found, declare ambiguous randomer (even if it has 'N's) as
    # non-ambiguous.
    for _, amb_bc in sorted(ambig_bcs):
        matches = False
        # Get list of accepted barcodes (sorted by decreasing frequency).
        # Sort is required each time as frequency changes when hits are
        # assigned from ambiguous to non-ambiguous randomer
        order_bcs = sorted([(len(hits), bc) for bc, hits in by_bc.items() if
                            bc in nonambig_bcs], reverse=True)

        for _, barcode in order_bcs:
            if _match(barcode, amb_bc, mismatches):
                matches = True
                by_bc[barcode].extend(by_bc.pop(amb_bc))
                break
        if not matches:
            nonambig_bcs.add(amb_bc)

    # Step #2: identify and accept randomers with strong support
    # Randomers that are supported by a threshold number of hits are accepted as unique.
    # Threshold is defined as the proportion (ratio_th) of the number of reads assigned
    # to the most frequent randomer.
    order_bcs = [(len(hits), bc) for bc, hits in by_bc.items()]
    order_bcs = sorted(order_bcs, reverse=True)
    min_hit_count = max(1, math.floor(order_bcs[0][0] * ratio_th))

    # Step #3: merge remaining
    # For each barcode with number of reads below the threshold, identify if there
    # exists any similar one. If there is, join the hits form second barcode to the
    # first one.
    merged = True
    while merged:
        # start with most frequent randomers first
        order_bcs = [(len(hits), bc) for bc, hits in by_bc.items()]
        order_bcs = sorted(order_bcs, reverse=True)
        merged = False
        for i, (_, barcode) in enumerate(order_bcs):
            for nhits2, barcode2 in order_bcs[i + 1:]:
                if nhits2 >= min_hit_count:
                    # do not try to merge barcodes with large support
                    continue
                if _match(barcode, barcode2, mismatches):
                    merged = True
                    by_bc[barcode].extend(by_bc.pop(barcode2))
            if merged:
                break


def _collapse(xlink_pos, by_bc, group_by, multimax=1):
    """
    Report number of cDNAs and reads in cross-link site on xlink_pos.

    Input parameter `by_bc` has te following structure:
    by_bc = {
        'AAA': [(middle_pos, end_pos, read_len, num_mapped, cigar, second_start),  # hit1
                (middle_pos, end_pos, read_len, num_mapped, cigar, second_start),  # hit2
                (middle_pos, end_pos, read_len, num_mapped, cigar, second_start),  # ...
        ]
        'AAT': [(middle_pos, end_pos, read_len, num_mapped, cigar, second_start),  # hit1
                (middle_pos, end_pos, read_len, num_mapped, cigar, second_start),  # hit2
                ]

    Counting the number of reads is easy - just count the number of hits per
    cross-link site.

    Counting the number of cDNAs is also easy - just count the number of
    different barcodes. However, following scenarions also need to be handled:

        * one read ban be mapped to multiple sites. In this case, the
          "contribution" of such read has to be divided equally to all positions
          that it maps to.
        * Also, longer reads should have proportionally greater "contribution".
          than the short ones.

    Upper two scenarions imply that each read contributes::

        weight = 1 * 1/a * b/c
        # a = number of hits
        # b = read length
        # c = sum(read lengths per same barcode)

    Another factor to take into account is also the possibility that a group of
    reads with equal start position and barcode represents multiple cross-links.
    Imagine a read starting 10 bp before exon-intron junction. One group of
    reads maps in the intron section and other reads skip the intron and map on
    next exon with the second part of the read. This can be solved by grouping
    by "second_start", which is the coordinate of the first nucleotide of the
    second part of the read. Each group with unique second start is treated as
    an independent cross-link event. This is done in function
    ``_separate_by_second_starts``

    Returns an object ``counts``::

        counts = {
            position: [cDNA_count, reads_count],
            123: [3.14, 42],
            124: [5.79, 16],
            ...
        }

    Parameters
    ----------
    xlink_pos : int
        Cross link position (genomic coordinate).
    by_bc : dict
        Dict with hits for each barcode.
    group_by : str
        Report by start, middle or end position.
    multimax : int
        Ignore reads, mapped to more than ``multimax`` places.

    Returns
    -------
    dict
        Number of cDNA and reads for each position.

    """
    group_by_index = ['start', 'middle', 'end'].index(group_by)

    # Container for cDNA and read counts:
    counts = {}

    for hits in by_bc.values():

        # separate in groups by second-start
        ss_groups = {}
        for read in hits:
            ss_groups.setdefault(read[4], []).append(read)

        for ss_group in ss_groups.values():

            # Sum of all read lengths per ss_group:
            sum_len_per_barcode = sum([i[2] for i in ss_group if i[3] <= multimax])

            for middle_pos, end_pos, read_len, num_mapped, _ in ss_group:
                if num_mapped > multimax:
                    continue
                grp_pos = (xlink_pos, middle_pos, end_pos)[group_by_index]
                weight = read_len / (num_mapped * sum_len_per_barcode)

                current_values = counts.get(grp_pos, (0, 0))
                upadated_values = (current_values[0] + weight, current_values[1] + 1)
                counts[grp_pos] = upadated_values

    return counts


def _intersects_with_annotaton(second_start, segmentation, chrom, strand):
    """
    Test if second_start corresopnds to any entry in segmentation.

    Returns
    -------
        bool
    Does the read's second_start corresopnd to any known segment in segmentation

    """
    for gene_content in segmentation.values():
        for transcript_id, transcript_content in gene_content.items():
            if transcript_id == 'gene_segment':
                continue
            for segment in transcript_content:
                if strand == '+':
                    if second_start == segment.start:
                        return True
                else:
                    if second_start == segment.stop:
                        return True
    return False


def _second_start(read, poss, strand, chrom, segmentation, holesize_th):
    """
    Return the coordinate of second start.

    If read is not split or we wish algorithm
    to think of read as linear, second_start equals to 0.
    """
    holes = [j - i - 1 for i, j in zip(poss, poss[1:])]
    # Get the size of the biggest hole:
    biggest_hole_size = max(holes) if holes else 0

    second_start = 0
    is_strange = False
    if not segmentation:
        # Effectively this means, that read is considered as it has no holes.
        if biggest_hole_size > holesize_th:
            # Still, read is not treated as on with distinct second_start.
            # However, it is reported as starnge:
            is_strange = True
    else:
        biggest_hole_size_index = holes.index(biggest_hole_size)
        # Take right border of hole on "+" and left border on "-" strand:
        if strand == '+':
            second_start = poss[biggest_hole_size_index + 1]
        else:
            second_start = poss[biggest_hole_size_index]

        # Read is strange if:
        # it is not intersecting with segmentation AND
        # if there actually is a hole
        if not _intersects_with_annotaton(second_start, segmentation, chrom, strand) and \
                biggest_hole_size != 0:
            is_strange = True

    return second_start, is_strange


def _get_read_data(read, metrics, mapq_th, segmentation=None, gap_th=4):
    """Extract neccessary data from read."""
    # NH (number of reported alignments) tag is required:
    if not read.has_tag('NH'):
        raise ValueError('"NH" tag not set for record: {}'.format(read.query_name))
    num_mapped = read.get_tag('NH')

    # Extract randomer sequence (barcode) from querry name (= read name)
    barcode = _get_random_barcode(read.query_name, metrics)
    metrics.bc_cn[barcode] = metrics.bc_cn.get(barcode, 0) + 1

    # position of cross-link is one nucleotide before start of read
    poss = read.get_reference_positions()
    if read.is_reverse:
        strand = '-'
        xlink_pos = poss[-1] + 1
        end_pos = poss[0]
    else:
        strand = '+'
        xlink_pos = poss[0] - 1
        xlink_pos = 1 if xlink_pos < 1 else xlink_pos  # Case of neg. pos on circular MT
        end_pos = poss[-1]

    chrom = read.reference_name
    second_start, is_strange = _second_start(read, poss, strand, chrom, segmentation, gap_th)
    if is_strange:
        metrics.strange_recs += 1

    # Position of middle nucleotide. Because we can have spliced reads, middle position is
    # not necessarily the middle of region the read maps to. Take one nucleotide upstream
    # of center in case length is even, happens by default on + strand
    idx = read.query_length // 2 - 1 if (strand == '-' and read.query_length % 2 == 0) \
        else read.query_length // 2
    middle_pos = poss[idx]

    return (xlink_pos, barcode, is_strange, strand, middle_pos, end_pos, read.query_length,
            num_mapped, second_start)


def _processs_bam_file(bam_fname, metrics, mapq_th, skipped, segmentation=None, gap_th=4):
    """
    Extract data from BAM file into chunks of genome.

    Parameters
    ----------
    bam_fname : str
        BAM file with mapped reads.
    metrics : iCount.Metrics
        Metrics object for storing analysis metadata.
    mapq_th : int
        Ignore hits with MAPQ < mapq_th.
    skipped : str
        Output BAM file to store reads that do not map as expected by segmentation and
        reference genome sequence. If read's second start does not fall on any of
        segmentation borders, it is considered problematic. If segmentation is not provided,
        every read in two parts with gap longer than gap_th is not used (skipped).
        All such reads are reported to the user for further exploration.
    segmentation : str
        File with segmentation (obtained by ``iCount segment``).
    gap_th : int
        Reads with gaps less than gap_th are treated as if they have no gap.

    Returns
    -------
    dict
        Internal structure of BAM file, described in docstring.
    list
        BAM file with

    """
    metrics.all_recs = 0  # All records
    metrics.notmapped_recs = 0  # Not mapped records
    metrics.mapped_recs = 0  # Mapped records
    metrics.lowmapq_recs = 0  # Records with insufficient quality
    metrics.used_recs = 0  # Records used in analysis (all - unmapped - lowmapq)
    metrics.invalidrandomer_recs = 0  # Records with invalid randomer
    metrics.norandomer_recs = 0  # Records with no randomer
    metrics.bc_cn = {}  # Barcode counter
    metrics.strange_recs = 0  # Strange records (not expected by segmentation)

    def finalize(reads_pending_fwd, reads_pending_rev, start, chrom, progress):
        """Yield appropriate data."""
        reads_to_process_fwd = {}
        for pos in list(reads_pending_fwd):
            if pos < start:
                reads_to_process_fwd[pos] = reads_pending_fwd.pop(pos)
        if reads_to_process_fwd:
            yield ((chrom, '+'), progress, reads_to_process_fwd)

        reads_to_process_rev = {}
        for pos in list(reads_pending_rev):
            if pos < start:
                reads_to_process_rev[pos] = reads_pending_rev.pop(pos)
        if reads_to_process_rev:
            yield ((chrom, '-'), progress, reads_to_process_rev)

    # Ensure sorted and and indexed input BAM file:
    LOGGER.info('Ensuring that bam file is sorted and indexed...')
    tmp_file = get_temp_file_name()
    pysam.sort('-o', tmp_file, bam_fname)  # pylint: disable=no-member
    pysam.index(tmp_file)  # pylint: disable=no-member
    genome_done = 0
    ann_data = None
    LOGGER.info('Detecting cross-links...')
    with AlignmentFile(tmp_file, 'rb') as bamfile:
        strange_bam = AlignmentFile(skipped, 'wb', header=bamfile.header)
        genome_size = sum([contig['LN'] for contig in bamfile.header['SQ']])
        for chrom in bamfile.references:
            chrom_len = bamfile.header['SQ'][bamfile.get_tid(chrom)]['LN']
            if segmentation:
                # pylint: disable=protected-access
                ann_data = iCount.genomes.segment._prepare_segmentation(segmentation, chrom)

            reads_pending_fwd = {}
            reads_pending_rev = {}
            read = None
            for read in bamfile.fetch(chrom):
                metrics.all_recs += 1
                if read.is_unmapped:
                    metrics.notmapped_recs += 1
                    continue
                metrics.mapped_recs += 1
                if read.mapping_quality < mapq_th:
                    metrics.lowmapq_recs += 1
                    continue
                metrics.used_recs += 1

                rdata = _get_read_data(
                    read, metrics, mapq_th, segmentation=ann_data, gap_th=gap_th)
                (xlink_pos, barcode, is_strange, strand), read_data = rdata[0:4], rdata[4:]

                if is_strange:
                    strange_bam.write(read)
                else:
                    if strand == '+':
                        reads_pending_fwd.setdefault(
                            xlink_pos, {}).setdefault(barcode, []).append(read_data)
                    else:
                        reads_pending_rev.setdefault(
                            xlink_pos, {}).setdefault(barcode, []).append(read_data)

            # Sliding window start (smaller coordinate)
            start = 0 if read is None else (0 if not read.positions else read.positions[0])
            progress = round(min((genome_done + start) / genome_size, 1.0), 4)

            for data in finalize(reads_pending_fwd, reads_pending_rev, start, chrom, progress):
                yield data

            start = chrom_len
            progress = round(min((genome_done + start) / genome_size, 1.0), 4)
            for data in finalize(reads_pending_fwd, reads_pending_rev, start, chrom, progress):
                yield data

            genome_done += chrom_len

    # Clean up:
    os.remove(tmp_file)

    # Report:
    LOGGER.info('All records in BAM file: %d', metrics.all_recs)
    LOGGER.info('Reads not mapped: %d', metrics.notmapped_recs)
    LOGGER.info('Mapped reads records (hits): %d', metrics.mapped_recs)
    LOGGER.info('Hits ignored because of low MAPQ: %d', metrics.lowmapq_recs)
    LOGGER.info('Records used for quantification: %d', metrics.used_recs)
    LOGGER.info('Records with invalid randomer info in header: %d', metrics.invalidrandomer_recs)
    LOGGER.info('Records with no randomer info: %d', metrics.norandomer_recs)
    LOGGER.info('Ten most frequent randomers:')
    top10 = sorted(
        [(count, barcode) for barcode, count in metrics.bc_cn.items()], reverse=True)[:10]
    for count, barcode in top10:
        LOGGER.info('    %s: %d', barcode, count)
    LOGGER.info('There are %d reads with second-start not falling on segmentation. They are '
                'reported in file: %s', metrics.strange_recs, skipped)


def run(bam, sites_unique, sites_multi, skipped, group_by='start', quant='cDNA',
        segmentation=None, mismatches=1, mapq_th=0, multimax=50, gap_th=4, ratio_th=0.1,
        report_progress=False):
    """
    Identify and quantify cross-linked sites.

    Interpret mapped sites and generate BED file with coordinates and
    number of cross-linked events.

    MAPQ is calculated mapq=int(-10*log10(1-1/Nmap)). By default we set
    the mapq_th to 0 to include all reads. Mapq score is very useful,
    because values coming from STAR are from a very limited set: 0 (5 or
    more multiple hits), 1 (4 or 3 multiple hits), 3 (2 multiple hits),
    255 (unique hit)

    Parameters
    ----------
    bam : str
        Input BAM file with mapped reads.
    sites_unique : str
        Output BED6 file to store data from uniquely mapped reads.
    sites_multi : str
        Output BED6 file to store data from multi-mapped reads.
    skipped : str
        Output BAM file to store reads that do not map as expected by segmentation and
        reference genome sequence. If read's second start does not fall on any of
        segmentation borders, it is considered problematic. If segmentation is not provided,
        every read in two parts with gap longer than gap_th is not used (skipped).
        All such reads are reported to the user for further exploration.
    group_by : str
        Assign score of a read to either 'start', 'middle' or 'end' nucleotide.
    quant : str
        Report number of 'cDNA' or number of 'reads'.
    mismatches : int
        Reads on same position with random barcode differing less than
        ``mismatches`` are merged together, if their ratio is below ratio_th.
    segmentation : str
        File with custon segmentation format (obtained by ``iCount segment``).
    mapq_th : int
        Ignore hits with MAPQ < mapq_th.
    multimax : int
        Ignore reads, mapped to more than ``multimax`` places.
    report_progress : bool
        Switch to report progress.
    gap_th : int
        Reads with gaps less than gap_th are treated as if they have no gap.
    ratio_th : float
        Ratio between the number of reads supporting a randomer versus the
        number of reads supporting the most frequent randomer. All randomers
        above this threshold are accepted as unique. Remaining are merged
        with the rest, allowing for the specified number of mismatches.

    Returns
    -------
    iCount.Metrics
        Metrics object, storing analysis metadata.

    """
    iCount.log_inputs(LOGGER, level=logging.INFO)  # pylint: disable=protected-access

    assert sites_unique.endswith(('.bed', '.bed.gz'))
    assert sites_multi.endswith(('.bed', '.bed.gz'))
    assert skipped.endswith(('.bam'))
    assert quant in ['cDNA', 'reads']
    assert group_by in ['start', 'middle', 'end']

    metrics = iCount.Metrics()

    unique, multi = {}, {}
    progress = 0
    for (chrom, strand), new_progress, by_pos in _processs_bam_file(
            bam, metrics, mapq_th, skipped, segmentation, gap_th):
        if report_progress:
            # pylint: disable=protected-access
            progress = iCount._log_progress(new_progress, progress, LOGGER)

        unique_by_pos = {}
        multi_by_pos = {}
        for xlink_pos, by_bc in by_pos.items():

            _merge_similar_randomers(by_bc, mismatches, ratio_th=ratio_th)

            # count uniquely mapped reads only
            _update(unique_by_pos, _collapse(xlink_pos, by_bc, group_by, multimax=1))
            # count all reads mapped les than multimax times
            _update(multi_by_pos, _collapse(xlink_pos, by_bc, group_by, multimax=multimax))

        unique.setdefault((chrom, strand), {}).update(unique_by_pos)
        multi.setdefault((chrom, strand), {}).update(multi_by_pos)

    # Write output
    val_index = ['cDNA', 'reads'].index(quant)
    _save_dict(unique, sites_unique, val_index=val_index)
    LOGGER.info('Saved to BED file (uniquely mapped reads): %s', sites_unique)
    _save_dict(multi, sites_multi, val_index=val_index)
    LOGGER.info('Saved to BED file (multi-mapped reads): %s', sites_multi)

    return metrics
