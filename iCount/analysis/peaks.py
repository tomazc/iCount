"""
Peak finding
------------

Determine local clusters of significantly cross-linked sites.

Read annotation and bedGraph with cross-linked sites. Perform permutation
analysis to determine significance (FDR) of individual sites. Return bedGraph of
significant sites.
"""
import math
import bisect
import logging

import numpy

import iCount

from collections import Counter

from iCount.files.bed import _f2s

LOGGER = logging.getLogger(__name__)


def sum_within_window(pos_val, hw=3):
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
        i_start = bisect.bisect_left(poss, p - hw, lo=max(0, i - hw), hi=i)
        i_stop = bisect.bisect_left(poss, p + 1 + hw,
                                    lo=i_start, hi=min(i + hw + 1, max_i))
        ret_list[inds[i]] = (p, sum(vals[i_start:i_stop]))
    return ret_list


def sum_within_window_nopos(pos_val, hw=3):
    if not pos_val:
        return []
    pos_val = sorted(pos_val)
    poss, vals = zip(*pos_val)

    ret_list = []
    max_i = len(poss)
    for i, p in enumerate(poss):
        i_start = bisect.bisect_left(poss, p - hw, lo=max(0, i - hw), hi=i)
        i_stop = bisect.bisect_left(poss, p + 1 + hw,
                                    lo=i_start, hi=min(i + hw + 1, max_i))
        ret_list.append(sum(vals[i_start:i_stop]))
    return ret_list


def cumulative_prob(vals, max_val):
    freqs, _ = numpy.histogram(vals, bins=range(max_val+1+1), density=True)
    freqs_cum = numpy.cumsum(freqs[::-1])
    return freqs_cum[::-1]


def cum_prob_within_window(pos_val, total_hits, hw=3):
    # extend counts to neighboring +-w nucleotides
    vals_extended = sum_within_window(pos_val, hw=hw)
    # calculate cumulative probabilities
    return cumulative_prob([v for _, v in vals_extended], total_hits), vals_extended


def cum_prob_within_window_nopos(pos_val, total_hits, hw=3):
    # extend counts to neighboring +-w nucleotides
    vals_extended = sum_within_window_nopos(pos_val, hw=hw)

    # calculate cumulative probabilities
    return cumulative_prob(vals_extended, total_hits)


# def get_rnd_distrib(size, total_hits, hw, perms=100):
#     """The simplest (and fastest) permutation test.
#
#     Not used by default.
#     """
#     rnd_cns = numpy.zeros(total_hits+1)
#     for i in range(perms):
#         rnd_hits = Counter(numpy.random.randint(size, size=total_hits))
#         rnd_hits_extended = sum_within_window_nopos(rnd_hits.items(), hw)
#         for v in rnd_hits_extended:
#             rnd_cns[v] += 1
#
#     s = sum(rnd_cns)
#     if s > 0:
#         rnd_ps = rnd_cns / s
#     else:
#         rnd_ps = rnd_cns
#
#     cum_prob_ret = numpy.cumsum(rnd_ps[::-1])
#     return cum_prob_ret[::-1]


ps_cache = {}


def get_avg_rnd_distrib(size, total_hits, hw, perms=100):
    """Return background distribution of peak heights for given region size
    and number of hits.

    We follow the modified FDR for peak height, proposed by:

    Yeo, G.W. et al. An RNA code for the FOX2 splicing regulator revealed
    by mapping RNA-protein interactions in stem cells. Nat. Struct. Mol.
    Biol. 16, 130–137 (2009).
    """
    global ps_cache

    c_key = (size, total_hits, hw, perms)
    if c_key in ps_cache:
        return ps_cache[c_key]

    rnd_ps = numpy.zeros((perms, total_hits+1))
    for pi in range(perms):
        rnd_hits = Counter(numpy.random.randint(size, size=total_hits))
        rnd_ps[pi, :] = cum_prob_within_window_nopos(rnd_hits.items(),
                                                     total_hits, hw=hw)

    cum_prob_ret = numpy.mean(rnd_ps, axis=0) + numpy.std(rnd_ps, axis=0)
    ps_cache[c_key] = cum_prob_ret
    return cum_prob_ret


def run(annotation, sites, peaks, scores=None, hw=3, fdr=0.05, perms=100,
        rnd_seed=42, features=['gene'], report_progress=False, group_by='gene_id'):
    """
    Calculate FDR of interaction at each cross-linked site.

    Parameters
    ----------
    annotation : str
        Path to input GTF annotation file.
    sites : str
        Path to input BED6 file listing all cross-linked sites.
    peaks : str
        Path to output BED6 file listing significant sites.
    scores : str
        Path to output tab-delimited file with detailed info on significance calculations.
    hw : int
        Half-window size.
    fdr : float
        FDR threshold.
    perms : int
        Number of permutations.
    rnd_seed : int
        Seed for random generator.
    features : list_str
        Features over which to perform analysis.
    report_progress : bool
        Print analysis progress.

    Returns
    -------
    None
        None

    """
    iCount.log_inputs(LOGGER, level=logging.INFO)

    # load annotation
    annotation = iCount.files.gtf.load(annotation)
    sites = iCount.files.bed.load(sites).sort().saveas()

    # assign cross-linked sites to regions
    LOGGER.info('Calculating intersection between annotation and cross-link file...')
    overlaps = annotation.intersect(sites, sorted=True, s=True, wo=True).saveas()
    hits_by_name = {}
    name_sizes = {}
    features_skipped_cn = 0
    LOGGER.info('Processing overlaps')
    for feature in overlaps:
        if feature[2] not in features:  # check if of correct feature type
            features_skipped_cn += 1
            continue
        group_id = feature.attrs[group_by]
        chrom = feature.chrom
        start = feature.start
        end = feature.stop
        assert start < end
        name = feature.name
        strand = feature.strand
        site_pos = int(feature.fields[10])
        site_score = float(feature.fields[13])
        hits_by_name.setdefault((chrom, strand, group_id, name), []).append((site_pos, site_score))
        name_sizes.setdefault((chrom, strand, group_id, name), set()).add((start, end))

    numpy.random.seed(rnd_seed)
    # calculate total length of each region
    name_sizes = dict(
        ((n, sum([e-s for s, e in v])) for n, v in name_sizes.items())
    )
    # calculate and assign FDRs to each cross-linked site

    # key is chrome, item is dict,
    # where key is (pos, strand), item is list of regions and associated pvals
    out_recs_scores = {}
    hits_by_name = sorted(hits_by_name.items())
    all_recs = len(hits_by_name)
    progress = 0
    while hits_by_name:
        gid, hits = hits_by_name.pop()
        region_size = name_sizes[gid]
        chrom, strand, group_id, name = gid
        region_hits = math.ceil(sum([v for _, v in hits]))

        if report_progress:
            new_progress = 1.0 - len(hits_by_name) / all_recs
            progress = iCount._log_progress(new_progress, progress, LOGGER)

        true_ps, hits_extended = cum_prob_within_window(hits, region_hits, hw)
        rnd_ps = get_avg_rnd_distrib(region_size, region_hits, hw, perms=perms)
        assert len(rnd_ps) == len(rnd_ps)
        val2fdr = [min(1.0, r/t) for r, t in zip(rnd_ps, true_ps)]
        for (p, val), (p2, val_extended) in zip(hits, hits_extended):
            assert p == p2
            score = val
            fdr_score = val2fdr[round(val_extended)]
            out_recs_scores.setdefault(chrom, {}).\
                setdefault((p, strand), []).\
                append((fdr_score, name, group_id, score, val_extended))

    peaks = iCount.files.gz_open(peaks, 'wt')
    if scores is not None:
        scores = iCount.files.gz_open(scores, 'wt')
        scores.write('chrome\tposition\tstrand\tname\tgroup_id\tscore'
                          '\tscore_extended\tFDR\n')

    for chrom, by_pos in sorted(out_recs_scores.items()):
        for (p, strand), annot_list in sorted(by_pos.items()):
            annot_list = sorted(annot_list)
            if scores:
                # all records are recorded in the score file
                for (fdr_score, name, group_id, score, val_extended) in annot_list:
                    # output in BED6 format:
                    # chrom, start, end, name, score, strand
                    scores.write(
                        '{:s}\t{:d}\t{:s}\t{:s}\t{:s}\t{:s}\t{:s}\t'
                        '{:.6f}\n'.format(chrom, p, strand, name, group_id,
                                          _f2s(score, dec=6),
                                          _f2s(val_extended, dec=6), fdr_score
                                          )
                    )

            # report minimum fdr_score for each position in BED6
            min_fdr_score = annot_list[0][0]
            if min_fdr_score < fdr:
                min_fdr_recs = [r for r in annot_list if r[0] == min_fdr_score]
                _, s_name, s_group_by, s_score, s_val_extended = zip(*min_fdr_recs)
                assert len(set(s_score)) == 1
                assert len(set(s_val_extended)) == 1  # this may not be true
                #  at borders of annotated regions
                score = s_score[0]
                name = ','.join(s_name)
                group_by = ','.join(s_group_by)
                o_str = '{:s}\t{:d}\t{:d}\t{:s}\t{:s}\t{:s}\t{:s}'.format(
                    chrom, p, p + 1, name, group_by, _f2s(score), strand
                )
                peaks.write('{:s}\n'.format(o_str))

    peaks.close()
    if scores is not None:
        scores.close()

    LOGGER.info('Bed file with significant peaks saved to %s', peaks)
    if scores:
        LOGGER.info('Scores for each cross-linked position saved to %s', scores)

    LOGGER.info('Done.')
