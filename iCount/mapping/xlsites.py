"""
Identify and quantify cross-linked sites
----------------------------------------

Transforms information in BAM file on mapped reads and their randomers into BED file.


Reads information on mapped reads and their associated randomer sequence
to produce BED file with cross-link sites.

Reads are first grouped by the start position. Multiple reads that start on
same position and have the same (similar) randomer are collapsed into one (
most common, if many of same frequency, take longest) read.

Collapsed reads can then be grouped in three different ways: by the start,
middle or end position.

When grouping by the start, we report one positions before read start. That
is the most likely cross-linked site.

Grouping by middle and end positions can be used for diagnostic purposes.


TODO: check overlap between unique and multimap BED files, should be small,
otherwise, we should think of a more approapriate quantification of (division
of randomers among) unique mapped sites that overlap with multimapped reads

"""

import pysam
import logging
import iCount


LOGGER = logging.getLogger(__name__)


def _match(s1, s2, mismatches):
    """
    Do sequence `s1` and `s2` have less or equal than ``mismatches``?

    Parameters
    ----------
    s1 : str
        First sequence.
    s2 : str
        Second sequence.
    mismatches : int
        Number of allowed mismatches between given sequences.

    Returns
    -------
    bool
        Do sequence `s1` and `s2` have less or equal than ``mismatches``

    """
    s1, s2 = s1.upper(), s2.upper()
    cn = sum([(c1 == 'N' or c2 == 'N' or c1 == c2) for c1, c2 in zip(s1, s2)])
    return max(len(s1), len(s2)) - cn <= mismatches


def _update(cur_vals, to_add):
    """
    Add the values from `to_add` to appropriate place in `cur_vals`

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
        prev_vals = cur_vals.get(pos, [0]*len(vals_to_add))
        cur_vals[pos] = [p + n for p, n in zip(prev_vals, vals_to_add)]


def _merge_similar_randomers(by_bc, mismatches):
    """
    Merge randomers on same site that are max ``mismatches`` different

    Input parameter `by_bc` has te following structure:
    by_bc = {
        'AAA': [(middle_pos, end_pos, read_len, num_mapped),  # hit1
                (middle_pos, end_pos, read_len, num_mapped),  # hit2
                (middle_pos, end_pos, read_len, num_mapped),  # ...
        ]
        'AAT': [(middle_pos, end_pos, read_len, num_mapped),  # hit1
                (middle_pos, end_pos, read_len, num_mapped),  # hit2
        ]

    Steps in function:
        0. Indentify ambigious randomers ('N' characters in barcode)

        1. For each ambigious randomer, identify similar non-ambigious one. If
        match is found, move hits form ambigious to the non-ambigious one. If no
        match is found, declare ambigious randomer (one that has 'N's) as
        non-ambigious anyway.

        2. For each barcode, identify if there exists any similar one. If
        there is, join the hits form second barcode to the first one.


    TODO: Code should be improved in step #1. Instead of finding any match,
    match with least difference should be found. Check also the skipped unit
    test in tests/test_xlsites.py

    Parameters
    ----------
    by_bc : dict
        Dictionary of barcodes and their hits
    mismatches : int
        Reads on same position with random barcode differing less than
        ``mismatches`` are grouped together.

    Returns
    -------
    None
        None, since input `by_bc` is modified in-place

    """
    # assign ambigious randomers to unambigious randomers
    accepted_bcs = set()  # accepted_barcodes
    ambig_bcs = []  # ambigious_barcodes
    for barcode in by_bc.keys():
        Ns = barcode.count('N')
        if Ns == 0:
            accepted_bcs.add(barcode)
        else:
            ambig_bcs.append((Ns, barcode))

    # Step #1:
    # For each ambigious randomer, identify similar non-ambigious one. If
    # match is found, move hits form ambigious to the non-ambigious one. If no
    # match is found, declare ambigious randomer (even if it has 'N's) as
    # non-ambigious anyway.
    for _, amb_bc in sorted(ambig_bcs):
        matches = False
        # Get list of accepted barcodes (sorted by decreasing frequency).
        # Sort is required each time as frequency changes when hits are
        # assigned from ambigious to unambigious randomer
        order_bcs = sorted([(len(hits), bc) for bc, hits in by_bc.items() if
                            bc in accepted_bcs], reverse=True)

        for _, bc in order_bcs:
            if _match(bc, amb_bc, mismatches):
                matches = True
                by_bc[bc].extend(by_bc.pop(amb_bc))
                break
        if not matches:
            accepted_bcs.add(amb_bc)

    # Step #2:
    # For each barcode, identify if there exists any similar one. If
    # there is, join the hits form second barcode to the first one.
    merged = True
    while merged:
        # start with most frequent randomers first
        order_bcs = [(len(hits), bc) for bc, hits in by_bc.items()]
        order_bcs = sorted(order_bcs, reverse=True)
        merged = False
        for i, (_, bc) in enumerate(order_bcs):
            for _, bc2 in order_bcs[i + 1:]:
                if _match(bc, bc2, mismatches):
                    merged = True
                    by_bc[bc].extend(by_bc.pop(bc2))
            if merged:
                break


def _separate_by_second_starts(hits):
    """
    Separate reads from ``hits`` in groups with same second start.

    Elemennts of list ``hits`` are tuples with the following entries::

        (middle_pos, end_pos, read_len, num_mapped, cigar, second_start)

    Parameters
    ----------
    hits : list
        List of tuples containing read information (described in upper section)

    Returns
    -------
    dict
        Entries from hits separated in groups with unique second start value.
    """
    # Sort reads from ``hits``:
    #     * first by their second start (sixth column)
    #     * than by read size in decreasing order (third column)
    hits = sorted(hits, key=lambda x: (x[5], -x[2]))

    # Determine the number of second start groups:
    second_starts = set([read[5] for read in hits])

    # separate hits into second_start groups:
    second_start_groups = {}
    for read in hits:
        second_start_groups.setdefault(read[5], []).append(read)

    return second_start_groups


def _collapse(xlink_pos, by_bc, group_by, multimax=1):
    """
    Report number of cDNAs and reads in cross-link site on xlink_pos

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
    gi = ['start', 'middle', 'end'].index(group_by)

    # Container for cDNA and read counts:
    counts = {}

    for bc, hits in by_bc.items():

        ss_groups = _separate_by_second_starts(hits)
        for ss_group in ss_groups.values():

            # Sum of all read lengths per ss_group:
            sum_len_per_barcode = sum([i[2] for i in ss_group if i[3] <= multimax])

            for middle_pos, end_pos, read_len, num_mapped, _, _ in ss_group:
                if num_mapped > multimax:
                    continue
                grp_pos = (xlink_pos, middle_pos, end_pos)[gi]
                w = read_len / (num_mapped * sum_len_per_barcode)

                current_values = counts.get(grp_pos, (0, 0))
                upadated_values = (current_values[0] + w, current_values[1] + 1)
                counts[grp_pos] = upadated_values

    return counts


def _second_start(start_positon, cigar):
    """
    Get the coordinate of the first nucleotide in the second part of a read.

    Read thas is mapped in two parts will look something like so::

        |-----exon-----|-----intron-----|

          |----R2.1----|      |---R2.2---|
       ...01234567890123456789012345...

    The corresopnding CIGAR value is ``((0,13), (3,7), (0,11))``. If the
    coordinate of the first nucleotide in the first part is 1000, than the
    coordinate of the first nucleotide in the second part is 1020.
    """
    first_hole_index = next((i for i, (code, _) in enumerate(cigar) if code == 3), None)
    if first_hole_index is None:
        return 0
    else:
        return start_positon + sum([num_nucs for _, num_nucs in cigar[:first_hole_index]])


def _processs_bam_file(bam_fname, metrics, mapq_th):
    """
    Extract BAM file to a dictionary.

    The structire of dictionary is the following::

        grouped = {
            ('chr1', '+'): (
                xlink_pos1: {
                    (barcode1, [
                        (middle_pos, end_pos, read_len, num_mapped, cigar, second_start)
                        ...
                    ]),
                    (barcode2, [
                        (middle_pos, end_pos, read_len, num_mapped, cigar, second_start)
                        ...
                    ]),
                    ....
                }),
            ),
        }

    Parameters
    ----------
    bam_fname : str
        BAM file with mapped reads.
    mapq_th : int
        Ignore hits with MAPQ < mapq_th.
    metrics : iCount.Metrics
        Metrics object, storing analysis metadata.
    cigar : bool
        Wheather to include cigar values for each read or not.

    Returns
    -------
    dict
        Internal structure of BAM file, described in docstring.
    """
    try:
        bamfile = pysam.AlignmentFile(bam_fname, 'rb')
    except OSError as e:
        raise ValueError('Error opening BAM file: {:s}'.format(bam_fname))

    # sanity check
    assert all(bamfile.getrname(i) == rname
               for i, rname in enumerate(bamfile.references))

    _cache_bcs = {}

    # counters
    metrics.all_recs = 0
    metrics.all_recs = 0  # All records
    metrics.notmapped_recs = 0  # Not mapped records
    metrics.mapped_recs = 0  # Mapped records
    metrics.lowmapq_recs = 0  # Records with insufficient quality
    metrics.used_recs = 0  # Records used in analysis (all - unmapped - lowmapq)
    metrics.invalidrandomer_recs = 0  # Records with invalid randomer
    metrics.norandomer_recs = 0  # Records with no randomer
    metrics.bc_cn = {}  # Barcode counter

    # group by start
    grouped = {}
    valid_nucs = set('ATCGN')
    for r in bamfile:
        metrics.all_recs += 1
        if r.is_unmapped:
            metrics.notmapped_recs += 1
            continue

        metrics.mapped_recs += 1

        if r.mapq < mapq_th:
            metrics.lowmapq_recs += 1
            continue

        metrics.used_recs += 1

        # NH (number of reported alignments) tag is required:
        if r.has_tag('NH'):
            num_mapped = r.get_tag('NH')
        else:
            raise ValueError('"NH" tag not set for record: {}'.format(r.qname))

        # Extract randomer sequence - it should be part of read id:
        if ':rbc:' in r.qname:
            bc = r.qname.rsplit(':rbc:', 1)[1].split(':')[0]
        elif ':' in r.qname:
            bc = r.qname.rsplit(':', 1)[1]
            if set(bc) - valid_nucs:
                # invalid barcode characters
                bc = ''
                metrics.invalidrandomer_recs += 1
        else:
            bc = ''
            metrics.norandomer_recs += 1

        bc = _cache_bcs.setdefault(bc, bc)  # reduce memory consumption
        metrics.bc_cn[bc] = metrics.bc_cn.get(bc, 0) + 1

        # position of cross-link is one nucleotide before start of read
        poss = sorted(r.positions)
        if r.is_reverse:
            strand = '-'
            xlink_pos = poss[-1] + 1
            end_pos = poss[0]
        else:
            strand = '+'
            xlink_pos = poss[0] - 1
            end_pos = poss[-1]
        chrome = bamfile.references[r.tid]

        # middle position is position of middle nucleotide. Because we can have
        # spliced reads, middle position is not necessarily the middle of
        # region the read maps to.
        i = len(poss)
        if i % 2 == 0:
            i = i // 2
            if strand == '-':
                # take one nucleotide upstream of center, happens by default
                # on + strand
                i = i - 1
        else:
            i = i // 2
        middle_pos = poss[i]
        read_len = len(r.seq)

        read_data = (middle_pos, end_pos, read_len, num_mapped, r.cigar,
                     _second_start(xlink_pos, r.cigar))

        # store hit data in a "dict of dicts" structure:
        grouped.setdefault((chrome, strand), {}).\
            setdefault(xlink_pos, {}).\
            setdefault(bc, []).\
            append(read_data)
    bamfile.close()
    return grouped


def run(bam_fname, unique_fname, multi_fname, group_by='start', quant='cDNA',
        mismatches=2, mapq_th=0, multimax=50):
    """
    Interpret mapped sites and generate BED file with coordinates and
    number of cross-linked events.

    MAPQ is calculated mapq=int(-10*log10(1-1/Nmap)). By default we set
    the mapq_th to 0 to include all reads. Mapq score is very useful,
    because values coming from STAR are from a very limited set: 0 (5 or
    more multiple hits), 1 (4 or 3 multiple hits), 3 (2 multiple hits),
    255 (unique hit)

    Parameters
    ----------
    bam_fname : str
        BAM file with mapped reads.
    unique_fname : str
        File to store data from uniquely mapped reads
    multi_fname : str
        File to store data from multi-mapped reads
    group_by : str
        Group reads together by 'start', 'middle' or 'end' nucleotide.
    quant : str
        Report number of 'cDNA' or number of 'reads'.
    mismatches : int
        Reads on same position with random barcode differing less than
        ``mismatches`` are grouped together.
    mapq_th : int
        Ignore hits with MAPQ < mapq_th.
    multimax : int
        Ignore reads, mapped to more than ``multimax`` places.

    Returns
    -------
    iCount.Metrics
        Metrics object, storing analysis metadata.
    """
    iCount.log_inputs(LOGGER, level=logging.INFO)

    assert quant in ['cDNA', 'reads']
    assert group_by in ['start', 'middle', 'end']

    metrics = iCount.Metrics()
    grouped = _processs_bam_file(bam_fname, metrics, mapq_th)

    # collapse duplicates
    unique = {}
    multi = {}
    while grouped:
        (chrome, strand), by_pos = grouped.popitem()
        unique_by_pos = {}
        multi_by_pos = {}
        while by_pos:
            xlink_pos, by_bc = by_pos.popitem()

            _merge_similar_randomers(by_bc, mismatches)
            # count uniquely mapped reads only
            _update(unique_by_pos, _collapse(xlink_pos, by_bc, group_by,
                                             multimax=1))
            # count also multi mapped reads
            _update(multi_by_pos, _collapse(xlink_pos, by_bc, group_by,
                                            multimax=multimax))

        unique[(chrome, strand)] = unique_by_pos
        multi[(chrome, strand)] = multi_by_pos

    # generate BED with cross-linked positions
    val_index = ['cDNA', 'reads'].index(quant)

    LOGGER.info('All records in BAM file: %d', metrics.all_recs)
    LOGGER.info('Reads not mapped: %d', metrics.notmapped_recs)
    LOGGER.info('Mapped reads records (hits): %d', metrics.mapped_recs)
    LOGGER.info('Hits ignored because of low MAPQ: %d', metrics.lowmapq_recs)
    LOGGER.info('Records used for quantification: %d', metrics.used_recs)
    LOGGER.info('Records with invalid randomer info in header: %d', metrics.invalidrandomer_recs)
    LOGGER.info('Records with no randomer info: %d', metrics.norandomer_recs)
    LOGGER.info('Ten most frequent randomers:')
    top10 = sorted([(cn, bc) for bc, cn in metrics.bc_cn.items()], reverse=True)[:10]
    for cn, bc in top10:
        LOGGER.info('    %s: %d', bc, cn)

    iCount.files.bed.save_dict(unique, unique_fname, val_index=val_index)
    LOGGER.info('Saved to BED file (uniquely mapped reads): %s', unique_fname)
    iCount.files.bed.save_dict(multi, multi_fname, val_index=val_index)
    LOGGER.info('Saved to BED file (multi-mapped reads): %s', multi_fname)

    return metrics
