"""
Generate BED file from BAM.
---------------------------

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

import iCount

# description and parameters needed for the analysis
analysis_name = 'xlsites'
analysis_description_short = 'identify and quantify cross-linked sites'
analysis_description = 'Transforms information in BAM file on mapped reads ' \
                       'and their randomers into BED file.'

params_opt = [
    (
        'groupby', 'choice-single', ('start', ['start', 'middle', 'end']),
        False,
        'Cross-link is defined by either the start, end or middle '
        'position of mapped reads.'
    ),
    (
        'quant', 'choice-single', ('cDNA', ['cDNA', 'reads']),
        False,
        'Report number of unique cDNAs or number of mapped reads.'
    ),
    (
        'mismatches', 'int_range', (1, 0, 5), False,
        'Number of tolerated mismatches when comparing two randomer barcodes.'
    ),
    (
        'multimax', 'int_range', (50, 1, 300), False,
        'Consider only reads with this number of hits or fewer.'
    ),
    (
        'mapqth', 'int_range', (0, 0, 255), False,
        'Consider only reads with at least this mapping quality.'
    ),
]

params_pos = [
    (
        'bam', 'BAM', 'in', 1,
        'Data on mapped iCLIP reads.'
    ),
    (
        'unique', 'folder', 'out', 1,
        'Name of BED file to store data on identified cross-links sites '
        'obtained from uniquely mapped reads.'
    ),
    (
        'multi', 'folder', 'out', 1,
        'Name of BED file to store data on identified cross-links sites '
        'obtained from multi-mapped reads.'
    ),
]


def _match(s1, s2, allowed_mismatches):
    """
    Do sequence `s1` and `s2` have less or equal than `allowed_mismatches`?

    :param str s1: First sequence
    :param str s2: Second sequence
    :param int allowed_mismatches: Number of allowed mismatches
    :return:
    :rtype: bool
    """
    s1, s2 = s1.upper(), s2.upper()
    cn = sum([(c1 == 'N' or c2 == 'N' or c1 == c2) for c1, c2 in zip(s1, s2)])
    return max(len(s1), len(s2)) - cn <= allowed_mismatches


def _update(cur_vals, to_add):
    """
    Add the values from `to_add` to appropriate place in `cur_vals`

    Note: cur_vals is updated in place!

    :param dict cur_vals: dict to be updated
    :param dict to_add: dict with update data
    :return: None
    :rtype: None
    """
    for pos, vals_to_add in to_add.items():
        prev_vals = cur_vals.get(pos, [0]*len(vals_to_add))
        cur_vals[pos] = [p + n for p, n in zip(prev_vals, vals_to_add)]


def _merge_similar_randomers(by_bc, randomer_mismatches):
    """
    Merge randomers on same site that are max `randomer_mismatches` different

    Input parameter `by_bc` has te following structure:
    by_bc = {
        'AAA': [(middle_pos, end_pos, read_len, num_mapped),  # hit1
                (middle_pos, end_pos, read_len, num_mapped),  # hit2
                (middle_pos, end_pos, read_len, num_mapped)]  # ...
        'AAT': [(middle_pos, end_pos, read_len, num_mapped),  # hit1
                (middle_pos, end_pos, read_len, num_mapped)]  # hit2

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

    :param dict by_bc: dictionary of barcodes and their hits
    :param int randomer_mismatches: number of allowed mismatches
    :return: None, since input `by_bc` is modified in-place
    :rtype: None
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
            if _match(bc, amb_bc, randomer_mismatches):
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
                if _match(bc, bc2, randomer_mismatches):
                    merged = True
                    by_bc[bc].extend(by_bc.pop(bc2))
            if merged:
                break


def _collapse(xlink_pos, by_bc, report_by, multimax=1):
    """
    Report number of cDNAs and reads in cross-link site on xlink_pos

    Input parameter `by_bc` has te following structure:
    by_bc = {
        'AAA': [(middle_pos, end_pos, read_len, num_mapped),  # hit1
                (middle_pos, end_pos, read_len, num_mapped),  # hit2
                (middle_pos, end_pos, read_len, num_mapped)]  # ...
        'AAT': [(middle_pos, end_pos, read_len, num_mapped),  # hit1
                (middle_pos, end_pos, read_len, num_mapped)]  # hit2

    Counting the number of reads is easy - just count the number of hits per
    cross-link site.

    Counting the number of cDNAs is also easy - just count the number of
    different barcodes. However, this is not the case, when one read has been
    mapped to many sites. In this case, the "contribution" of such read has to
    be divided equally to all positions that it maps to. Also, longer reads,
    should have greater influence. This is why the weight for cDNA count is:

        weight = a / (b * c)

        a = read length
        b = sum(read-lengths per barcode position)
        c = number of hits

    :param int xlink_pos: cross link position (genomic coordinate)
    :param dict by_bc: dict with hits for each barcode
    :param str report_by: report by start, middle or end position
    :param int multimax: consider only reads with multimax hits or fewer
    """
    gi = ['start', 'middle', 'end'].index(report_by)

    cDNAs = {}
    reads = {}
    for bc, hits in by_bc.items():
        cdna_cn = {}  # cDNA count
        read_cn = {}  # read count
        # Sum of all read lengths per barcode:
        sum_len_per_barcode = sum([i[2] for i in hits if i[3] <= multimax])
        for middle_pos, end_pos, read_len, num_mapped in hits:
            if num_mapped > multimax:
                continue
            grp_pos = (xlink_pos, middle_pos, end_pos)[gi]
            w = read_len / (num_mapped * sum_len_per_barcode)
            cdna_cn[grp_pos] = cdna_cn.get(grp_pos, 0) + w
            read_cn[grp_pos] = read_cn.get(grp_pos, 0) + 1

        # update cDNAs supporting the site (by which we group by)
        for grp_pos, w in cdna_cn.items():
            cDNAs[grp_pos] = cDNAs.get(grp_pos, 0) + w

        # update number of reads supporting the site (by which we group by)
        for grp_pos, cn in read_cn.items():
            reads[grp_pos] = reads.get(grp_pos, 0) + cn

    # merge cDNAs and reads into a tuple
    retd = {}
    for grp_pos, tot_cDNA in cDNAs.items():
        tot_reads = reads[grp_pos]
        retd[grp_pos] = (tot_cDNA, tot_reads)
    return retd


def run(bam_fname, unique_fname, multi_fname, group_by='start', quant='cDNA',
        randomer_mismatches=2, mapq_th=0, multimax=50):
    """
    Interpret mapped sites and generate BED file with coordinates and
    number of cross-linked events.

    MAPQ is calculated mapq=int(-10*log10(1-1/Nmap)). By default we set
    the mapq_th to 0 to include all reads. Mapq score is very useful,
    because values coming from STAR are from a very limited set: 0 (5 or
    more multiple hits), 1 (4 or 3 multiple hits), 3 (2 multiple hits),
    255 (unique hit)

    :param str bam_fname: path to bam filename
    :param str unique_fname: file to store data from uniquely mapped reads
    :param str multi_fname: file to store data from multi-mapped reads
    :param str group_by:
    :param str quant:
    :param int randomer_mismatches:
    :param int mapq_th: Ignore hits with MAPQ lower than this threshold value.
    :param int multimax:
    :return:
    :rtype:
    """
    assert quant in ['cDNA', 'reads']
    assert group_by in ['start', 'middle', 'end']
    try:
        bamfile = pysam.AlignmentFile(bam_fname, 'rb')
    except OSError as e:
        raise ValueError('Error opening BAM file: {:s}'.format(bam_fname))

    # sanity check
    assert all(bamfile.getrname(i) == rname \
               for i, rname in enumerate(bamfile.references))

    # counters
    all_recs = 0  # All records
    notmapped_recs = 0  # Not mapped records
    mapped_recs = 0  # Mapped records
    lowmapq_recs = 0  # Records with insufficient quality
    used_recs = 0  # Records used in analysis (all - unmapped - lowmapq)
    invalidrandomer_recs = 0  # Records with invalid randomer
    norandomer_recs = 0  # Records with no randomer
    bc_cn = {}  # Barcode counter
    _cache_bcs = {}

    # group by start
    grouped = {}
    valid_nucs = set('ATCGN')
    for r in bamfile:
        all_recs += 1
        if r.is_unmapped:
            notmapped_recs += 1
            continue

        mapped_recs += 1

        if r.mapq < mapq_th:
            lowmapq_recs += 1
            continue

        used_recs += 1

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
                invalidrandomer_recs += 1
        else:
            bc = ''
            norandomer_recs += 1

        bc = _cache_bcs.setdefault(bc, bc)  # reduce memory consumption
        bc_cn[bc] = bc_cn.get(bc, 0) + 1

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

        # store hit data in a "dict of dicts" structure:
        grouped.setdefault((chrome, strand), {}).\
                setdefault(xlink_pos, {}).\
                setdefault(bc, []).\
                append((middle_pos, end_pos, read_len, num_mapped))
    bamfile.close()

    # collapse duplicates
    unique = {}
    multi = {}
    while grouped:
        (chrome, strand), by_pos = grouped.popitem()
        unique_by_pos = {}
        multi_by_pos = {}
        while by_pos:
            xlink_pos, by_bc = by_pos.popitem()

            _merge_similar_randomers(by_bc, randomer_mismatches)
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
    iCount.files.bed.save_dict(unique, unique_fname, val_index=val_index)
    iCount.files.bed.save_dict(multi, multi_fname, val_index=val_index)

    return all_recs, notmapped_recs, mapped_recs, lowmapq_recs, \
           used_recs, invalidrandomer_recs, norandomer_recs, bc_cn
