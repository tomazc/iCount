"""
Peak finding
------------

Read annotation and bedGraph with cross-linked sites.
Perform permutation analysis to determine significance (FDR) of individual
sites.
Return bedGraph of significant sites.

"""
from collections import Counter
import bisect
import numpy

import iCount

# description and parameters needed for the analysis
analysis_name = 'peaks'
analysis_description_short = 'peak analysis'
analysis_description = 'Determine local clusters of significantly ' \
                       'cross-linked sites by performing permutation ' \
                       'analysis.'

params_opt = [
    (
        'hw', 'int_range', (3, 0, 200), False,
        'Half-window size.'
    ),
    (
        'fdr', 'float_range', (0.05, 0.01, 0.25), False,
        'Significance threshold.'
    ),
    (
        'perms', 'int_range', (100, 10, 1000), False,
        'Number of random permutations needed to determine statistical '
        'significance.'
    ),
# TODO: ignore for now, until change of annotation format from BED6 to GTF.
#    (
#        'regions', 'str_list', ['gene'],  False,
#        'Calculate enrichment for cross-links within these types (as given
#    # in '
#        '3rd column in GTF).'
#    )
]

params_pos = [
    (
        'annotation', 'BED6', 'in', 1,
        '(input) BED6 file with gene regions.'
    ),
    (
        'sites', 'BED6', 'in', 1,
        '(input) BED6 file with cross-linked sites.'
    ),
    (
        'peaks', 'BED6', 'out', 1,
        '(output) BED6 with significant cross-linked sites.'
    ),
    (
        'scores', 'tab', 'out', 1,
        '(output) tab-separated table with scores for each cross-linked site.'
    ),
]


def sum_within_window(pos_val, w=3):
    """Sum counts within window placed on top of each site

    The returned list preserves the order of positions given on input.
    """
    if not pos_val:
        return []
    pos_val_ind = sorted((p, v, i) for i, (p, v) in enumerate(pos_val))
    poss, vals, inds = zip(*pos_val_ind)

    ret_list = [None]*len(pos_val)
    max_i = len(poss)
    for i, p in enumerate(poss):
        i_start = bisect.bisect_left(poss, p-w, lo=max(0, i-w), hi=i)
        i_stop = bisect.bisect_left(poss, p+1+w,
                                    lo=i_start, hi=min(i+w+1, max_i))
        ret_list[inds[i]] = (p, sum(vals[i_start:i_stop]))
    return ret_list


def sum_within_window_nopos(pos_val, w=3):
    if not pos_val:
        return []
    pos_val = sorted(pos_val)
    poss, vals = zip(*pos_val)

    ret_list = []
    max_i = len(poss)
    for i, p in enumerate(poss):
        i_start = bisect.bisect_left(poss, p-w, lo=max(0, i-w), hi=i)
        i_stop = bisect.bisect_left(poss, p+1+w,
                                    lo=i_start, hi=min(i+w+1, max_i))
        ret_list.append(sum(vals[i_start:i_stop]))
    return ret_list


def cumulative_prob(vals, max_val):
    freqs, _ = numpy.histogram(vals, bins=range(max_val+1+1), density=True)
    freqs_cum = numpy.cumsum(freqs[::-1])
    return freqs_cum[::-1]


def cum_prob_within_window(pos_val, total_hits, w=3):
    # extend counts to neighboring +-w nucleotides
    vals_extended = sum_within_window(pos_val, w=w)
    # calculate cumulative probabilities
    return cumulative_prob([v for _, v in vals_extended], total_hits), vals_extended


def cum_prob_within_window_nopos(pos_val, total_hits, w=3):
    # extend counts to neighboring +-w nucleotides
    vals_extended = sum_within_window_nopos(pos_val, w=w)

    # calculate cumulative probabilities
    return cumulative_prob(vals_extended, total_hits)


def diff_ps(ps1, ps2):
    return sum(abs(p1-p2) for p1, p2 in zip(ps1, ps2))


def get_rnd_distrib(size, total_hits, w, perms=100):
    """The simplest (and fastest) permutation test."""
    rnd_cns = numpy.zeros(total_hits+1)
    for i in range(perms):
        rnd_hits = Counter(numpy.random.randint(size, size=total_hits))
        rnd_hits_extended = \
            iCount.analysis.peaks.sum_within_window_nopos(rnd_hits.items(), w)
        for v in rnd_hits_extended:
            rnd_cns[v] += 1

    s = sum(rnd_cns)
    if s > 0:
        rnd_ps = rnd_cns / s
    else:
        rnd_ps = rnd_cns

    cum_prob_ret = numpy.cumsum(rnd_ps[::-1])
    return cum_prob_ret[::-1]


ps_cache = {}

def get_avg_rnd_distrib(size, total_hits, w, perms=100):
    """Return background distribution of peak heights for given region size
    and number of hits.

    We follow the modified FDR for peak height, proposed by:

    Yeo, G.W. et al. An RNA code for the FOX2 splicing regulator revealed
    by mapping RNA-protein interactions in stem cells. Nat. Struct. Mol.
    Biol. 16, 130–137 (2009).
    """
    global ps_cache

    c_key = (size, total_hits, w, perms)
    if c_key in ps_cache:
        return ps_cache[c_key]

    rnd_ps = numpy.zeros((perms, total_hits+1))
    for pi in range(perms):
        rnd_hits = Counter(numpy.random.randint(size, size=total_hits))
        rnd_ps[pi, :] = cum_prob_within_window_nopos(rnd_hits.items(),
                                                     total_hits, w=w)

    cum_prob_ret = numpy.mean(rnd_ps, axis=0) + numpy.std(rnd_ps, axis=0)
    ps_cache[c_key] = cum_prob_ret
    return cum_prob_ret


def run(fin_annotation, fin_sites, fout_peaks, fout_scores=None, hw=3,
        fdr=0.05, perms=100, rnd_seed=42): #, regions=['gene']):
    """Calculate FDR of interaction at each cross-linked site.

    """

    # load annotation
    annotation = iCount.files.gtf.load(fin_annotation)
    sites = iCount.files.bed.load(fin_sites).sort().saveas()

    # assign cross-linked sites to regions
    overlaps = annotation.intersect(sites, sorted=True, s=True, wo=True).saveas()
    hits_by_name = {}
    name_sizes = {}
#    skip_cn = 0
    for feature in overlaps:
#        if feature.type not in regions:
#            skip_cn += 1
#            continue
        chrom = feature.chrom
        start = feature.start
        end = feature.stop
        assert start < end
        name = feature.name
        strand = feature.strand
        site_pos = int(feature.fields[7])
        site_score = int(feature.fields[10])
        hits_by_name.setdefault((chrom, strand, name), []).append((site_pos,
                                                                   site_score))
        name_sizes.setdefault((chrom, strand, name), set()).add((start, end))

    numpy.random.seed(rnd_seed)
    # calculate total length of each region
    name_sizes = dict(
        ((n, sum([e-s for s, e in v])) for n, v in name_sizes.items())
    )
    # calculate and assign FDRs to each cross-linked site
    out_recs_scores = []
    hits_by_name = sorted(hits_by_name.items())
    all_recs = len(hits_by_name)
    cur_perc = 0
    while hits_by_name:
        gid, hits = hits_by_name.pop()
        region_size = name_sizes[gid]
        chrom, strand, name = gid
        region_hits = sum([v for _, v in hits])
        new_perc = '\r{:.1f}%'.format(100.0*(1.0-len(hits_by_name)/all_recs))
        if new_perc != cur_perc:
            cur_perc = new_perc
            print(cur_perc, end="", flush=True),
        true_ps, hits_extended = cum_prob_within_window(hits, region_hits, hw)
        rnd_ps = get_avg_rnd_distrib(region_size, region_hits, hw, perms=perms)
        assert len(rnd_ps) == len(rnd_ps)
        val2fdr = [r/t for r, t in zip(rnd_ps, true_ps)]
        for (p, val), (p2, val_extended) in zip(hits, hits_extended):
            assert p == p2
            score = str(val)
            fdr_score = val2fdr[val_extended]
            out_recs_scores.append((chrom, p, name, score, strand,
                                    val_extended, fdr_score))
    out_recs_scores.sort()

    fout_peaks = iCount.files.gz_open(fout_peaks, 'wt')
    if fout_scores is not None:
        fout_scores = iCount.files.gz_open(fout_scores, 'wt')

    for (chrom, p, name, score, strand, val_extended, fdr_score) in \
            out_recs_scores:
        # output in BED6 format:
        # chrom, start, end, name, score, strand
        o_str = '{:s}\t{:d}\t{:d}\t{:s}\t{:s}\t{:s}'.format(chrom, p,
                                                            p + 1, name,
                                                            score, strand)
        if fdr_score < fdr:
            fout_peaks.write('{:s}\n'.format(o_str))

        if fout_scores is not None:
            fout_scores.write('{:s}\t{:s}\t{:f}\n'.format(o_str,
                                                          str(val_extended),
                                                          fdr_score))
    fout_peaks.close()
    if fout_scores is not None:
        fout_scores.close()
    print('\ndone')
