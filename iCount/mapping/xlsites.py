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


def match(s1, s2, allowed_mismatches):
    cn = sum([(c1 == 'N' or c2 == 'N' or c1 == c2) for c1, c2 in zip(s1, s2)])
    return len(s1) - cn <= allowed_mismatches


def _update(cur_vals, to_add):
    for pos, vals_to_add in to_add.items():
        prev_vals = cur_vals.get(pos, [0]*len(vals_to_add))
        cur_vals[pos] = [p+n for p, n in zip(prev_vals, vals_to_add)]


def _merge_similar_randomers(by_bc, randomer_mismatches):
    # assign ambigious randomers to unambigious randomers
    accepted_bcs = set()
    ambig_bcs = []
    for bc, hits in by_bc.items():
        Ns = bc.count('N')
        if Ns == 0:
            accepted_bcs.add(bc)
        else:
            ambig_bcs.append((Ns, bc))

    # assign ambigious (in increasing order of Ns in randomer)
    ambig_bcs.sort()
    for _, amb_bc in ambig_bcs:
        matches = False
        # sort unambigious by decreasing frequency, need to sort each time as
        # frequency changes when hits are assigned from ambigious to
        # unambigious randomer
        order_bcs = [(len(hits), bc) for bc, hits in by_bc.items()]
        order_bcs = sorted(order_bcs, reverse=True)
        for bc in order_bcs:
            if bc not in accepted_bcs:
                continue
            if match(amb_bc, bc, randomer_mismatches):
                matches = True
                amb_hits = by_bc.pop(amb_bc)
                by_bc[bc].extend(amb_hits)
                break
        if not matches:
            accepted_bcs.add(amb_bc)

    # merge similar radomers first
    merged = True
    while merged:
        # start with most frequent randomers first
        order_bcs = [(len(hits), bc) for bc, hits in by_bc.items()]
        order_bcs = sorted(order_bcs, reverse=True)
        merged = False
        for i, (_, bc) in enumerate(order_bcs):
            for bc2 in order_bcs[i+1:]:
                if match(bc, bc2, randomer_mismatches):
                    merged = True
                    amb_hits = by_bc.pop(amb_bc)
                    by_bc[bc].extend(amb_hits)
            if merged:
                break

    return by_bc


def _collapse(xlink_pos, by_bc, report_by, multimax=1):
    """"""
    gi = ['start', 'middle', 'end'].index(report_by)

    cDNAs = {}
    reads = {}
    for bc, hits in by_bc.items():
        w_e = {}  # subdivide among positions that we group by
        w_d = 0.0
        r_cn = {}
        for middle_pos, end_pos, read_len, num_mapped in hits:
            if num_mapped > multimax:
                continue
            grp_pos = (xlink_pos, middle_pos, end_pos)[gi]
            w_e[grp_pos] = w_e.get(grp_pos, 0) + read_len/num_mapped
            w_d += read_len
            r_cn[grp_pos] = r_cn.get(grp_pos, 0) + 1

        # update cDNAs supporting the site (by which we group by)
        for grp_pos, e in w_e.items():
            w = e / w_d
            cDNAs[grp_pos] = cDNAs.get(grp_pos, 0) + w

        # update number of reads supporting the site (by which we group by)
        for grp_pos, cn in r_cn.items():
            reads[grp_pos] = reads.get(grp_pos, 0) + cn

    # merge cDNAs and reads into a tuple
    retd = {}
    for grp_pos, tot_cDNA in cDNAs.items():
        tot_reads = reads[grp_pos]
        retd[grp_pos] = (tot_cDNA, tot_reads)
    return retd


def run(bam_fname, unique_fname, multi_fname, group_by='start', quant='cDNA',
        randomer_mismatches=2, mapq_th=0, multimax=50):
    """Interpret mapped sites and generate BED file with coordinates and
    number of cross-linked events.

    MAPQ is calculated mapq=int(-10*log10(1-1/Nmap)). By default we set
    the mapq_th to 0 to include all reads. Mapq score is very useful,
    because values coming from STAR are from a very limited set: 0 (5 or
    more multiple hits), 1 (4 or 3 multiple hits), 3 (2 multiple hits),
    255 (unique hit)

    :param bam_fname:
    :param unique_fname:
    :param multi_fname:
    :param group_by:
    :param randomer_mismatches:
    :param mapq_th: Ignore hits with MAPQ lower than this value.
    :param multimax:
    :return:
    """
    assert quant in ['cDNA', 'reads']
    assert group_by in ['start', 'middle', 'end']
    try:
        bamfile = pysam.Samfile(bam_fname, 'rb')
    except:
        print('Error opening BAM file: {:s}'.format(bam_fname))
        return

    # sanity check
    assert all(bamfile.getrname(i) == rname \
               for i, rname in enumerate(bamfile.references))

    # counters
    all_recs = 0
    notmapped_recs = 0
    mapped_recs = 0
    lowmapq_recs = 0
    used_recs = 0
    invalidrandomer_recs = 0
    norandomer_recs = 0
    bc_cn = {}
    _cache_bcs = {}

    # group by start
    grouped = {}
    valid_nucs = set('ATCGN')
    for r in bamfile:
        all_recs += 1
        if r.is_unmapped:
            notmapped_recs += 1
            continue

        if r.mapq < mapq_th:
            lowmapq_recs += 1
            continue

        mapped_recs += 1

        # record will be used
        used_recs += 1

        # get number of times mapped
        if r.has_tag('NH'):
            num_mapped = r.get_tag('NH')
        else:
            # we require each record to have the NH tag - number of reported
            # alignments
            print('Error, NH tags are required in BAM file.')
            return

        # randomer must be part of read id
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

        # store hit data
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
