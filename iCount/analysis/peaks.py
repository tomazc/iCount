"""
Peak finding
------------

Read annotation and bedGraph with cross-linked sites.
Perform permutation analysis to determine significance (FDR) of individual
sites.
Return bedGraph of significant sites.

"""
import bisect
import random
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
        'size', 'int_range', (3, 1, 200), False,
        'Peak width.'
    ),
    (
        'fdr', 'float_range', (0.05, 0.01, 0.25), False,
        'Significance threshold.'
    ),
    (
        'rnd', 'int_range', (100, 10, 250), False,
        'Number of random permutations needed to determine statistical '
        'significance.'
    ),
    (
        'regions', 'choice-single', ('one', ['one', 'separately']),  False,
        'Determine significance based on entire gene or withing each '
        'segment (exons, introns, UTRs, ...).'
    )
]

params_pos = [
    (
        'annotation', 'GTF', 'in',
        '(input) GTF file with gene models.'
    ),
    (
        'sites', 'bedGraph', 'in',
        '(input) bedGraph with cross-linked sites.'
    ),
    (
        'peaks', 'bedGraph', 'out',
        '(output) bedGraph with significant cross-linked sites.'
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
        i_start = bisect.bisect_left(poss, p-w, lo=max(0,i-w), hi=i)
        i_stop = bisect.bisect_left(poss, p+1+w, lo=i_start, hi=min(i+w+1,
                                                                    max_i))
        s = sum(vals[i_start:i_stop])
        ret_list[inds[i]] = (p, s)
    return ret_list


def sum_within_window_nopos(pos_val, w=3):
    if not pos_val:
        return []
    pos_val = sorted(pos_val)
    poss, vals = zip(*pos_val)

    ret_list = []
    max_i = len(poss)
    for i, p in enumerate(poss):
        i_start = bisect.bisect_left(poss, p-w, lo=max(0,i-w), hi=i)
        i_stop = bisect.bisect_left(poss, p+1+w, lo=i_start, hi=min(i+w+1,
                                                                    max_i))
        s = sum(vals[i_start:i_stop])
        ret_list.append(s)
    return ret_list


def cumulative_prob(vals):
    assert min(vals) >= 0
    max_val = max(vals)
    freqs, _ = numpy.histogram(vals, bins=range(max_val+1+1), density=True)
    freqs_cum = numpy.cumsum(freqs[::-1])
    freqs_cum = freqs_cum[::-1]
    return freqs_cum


def cum_prob_within_window(pos_val, w=3):
    # extend counts to neighboring +-w nucleotides
    vals_extended = sum_within_window(pos_val, w=w)
    # calculate cumulative probabilities
    return cumulative_prob([v for _, v in vals_extended]), vals_extended


def cum_prob_within_window_nopos(pos_val, w=3):
    # extend counts to neighboring +-w nucleotides
    vals_extended = sum_within_window_nopos(pos_val, w=w)

    # calculate cumulative probabilities
    return cumulative_prob(vals_extended)


# def get_rnd_distrib(size, total_hits, w, perms=100):
#     """The simplest (and fastest) permutation test."""
#     rnd_distrib = []
#     for _ in range(perms):
#         rnd_hits = {}
#         for _ in range(total_hits):
#             rnd_pos = random.randrange(size)
#             rnd_hits[rnd_pos] = rnd_hits.get(rnd_pos, 0) + 1
#
#         rnd_hits_extended = \
#             iCount.analysis.peaks.sum_within_window(list(rnd_hits.items()), w)
#
#         _, vals = zip(*rnd_hits_extended)
#         max_val = max(vals) + 1
#         missing_indexes = max_val - len(rnd_distrib)
#         if missing_indexes > 0:
#             rnd_distrib.extend([0]*missing_indexes)
#         for v in vals:
#             rnd_distrib[v] += 1
#
#     s = max(1, sum(rnd_distrib))
#     freqs = [cn/s for cn in rnd_distrib]
#     freqs_cum = numpy.cumsum(freqs[::-1])
#     freqs_cum = freqs_cum[::-1]
#     return list(freqs_cum)


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

    if size > 1000:
        c_size = size//1000
        c_size = c_size * 1000
    else:
        c_size = size
    if total_hits > 300:
        c_total_hits = total_hits // 50
        c_total_hits = c_total_hits * 50
    else:
        c_total_hits = total_hits // 10
        c_total_hits = c_total_hits * 10

    c_key = (c_size, c_total_hits, w, perms)
    print(size, total_hits, w, perms, c_key)
    if c_key in ps_cache:
        print('cache hit!')
        return ps_cache[c_key]

    rnd_distrib = []
    for _ in range(perms):
        rnd_hits = {}
        for _ in range(total_hits):
            rnd_pos = random.randrange(size)
            rnd_hits[rnd_pos] = rnd_hits.get(rnd_pos, 0) + 1

        pos_val = list(rnd_hits.items())
        ps = cum_prob_within_window_nopos(pos_val, w=w)

        missing_indexes = max(0, len(ps) - len(rnd_distrib))
        for i in range(missing_indexes):
            rnd_distrib.append([])

        for v, p in enumerate(ps):
            rnd_distrib[v].append(p)

    rnd_distrib_avg = []
    rnd_distrib_std = []
    cum_prob_ret = []
    for ps in rnd_distrib:
        if len(ps) < perms:
            ps.extend([0.0]*(perms-len(ps)))
        m = numpy.mean(ps)
        s = numpy.std(ps)
        rnd_distrib_avg.append(m)
        rnd_distrib_std.append(s)
        cum_prob_ret.append(m+s)

    ps_cache[c_key] = cum_prob_ret
    return cum_prob_ret


def run(fin_annotation, fin_sites, fout_peaks, fout_scores=None, hw=3,
        fdr=0.05, perms=100, regions='one'):
    """Calculate FDR of interaction at each cross-linked site.

    """

    # load annotation
    annotation = iCount.files.gtf.load(fin_annotation)
    sites = iCount.files.bed.load(fin_sites)

    # assign cross-linked sites to regions
    overlaps = annotation.intersect(sites, sorted=True, s=True, wo=True).saveas()
    hits_by_name = {}
    name_sizes = {}
    for feature in overlaps:
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

    if fout_peaks:
        fout_peaks = iCount.files.gz_open(fout_peaks, 'wt')
    if fout_scores:
        fout_scores = iCount.files.gz_open(fout_scores, 'wt')

    # calculate total length of each region
    name_sizes = dict(
        ((n, sum([e-s for s, e in v])) for n, v in name_sizes.items())
    )
    # calculate and assign FDRs to each cross-linked site
    for gid, hits in sorted(hits_by_name.items()):
        print(gid)
        size = name_sizes[gid]
        chrom, strand, name = gid
        total_hits = sum([v for _, v in hits])
        true_ps, hits_extended = cum_prob_within_window(hits, hw)
        rnd_ps = get_avg_rnd_distrib(size, total_hits, hw, perms=perms)
        rnd_ps.extend([0.0]*(len(true_ps)-len(rnd_ps)))
        val2fdr = [r/t for r, t in zip(rnd_ps, true_ps)]

        # output in BED6 format:
        # chrom, start, end, name, score, strand
        for (p, val), (p2, val_extended) in zip(hits, hits_extended):
            assert p == p2
            score = str(val)
            fdr_score = val2fdr[val_extended]
            o_str = '{:s}\t{:d}\t{:d}\t{:s}\t{:s}\t{:s}'.format(chrom, p,
                                                                p+1, name,
                                                                score, strand)
            if fdr_score < fdr:
                fout_peaks.write('{:s}\n'.format(o_str))
            if fout_scores:
                fout_scores.write('{:s}\t{:s}\t{:f}\n'.format(o_str,
                                                              str(val_extended),
                                                              fdr_score))

    fout_peaks.close()
    fout_scores.close()

#
# # code
# def extended_sum(hits, flank):
#     """ hits must be sorted by increasing position """
#     extended_hits = []
#     sum_from = 0
#     L = len(hits)
#     for pos, value in hits:
#         while sum_from < L and pos - hits[sum_from][0] > flank:
#             sum_from += 1
#         n_recs_min = pos
#         n_recs_max = pos
#         sum_val = 0.0
#         j = sum_from
#         while j < L and hits[j][0] - pos <= flank:
#             hit_pos, hit_val = hits[j]
#             sum_val += hit_val
#             if hit_pos < n_recs_min:
#                 n_recs_min = hit_pos
#             elif hit_pos > n_recs_max:
#                 n_recs_max = hit_pos
#             j += 1
#         extended_hits.append((pos, sum_val, n_recs_min, n_recs_max))
#     return extended_hits
#
#
# def cumulative_probs(pos_val_l):
#     val_freqs = {}
#     for pos, val, n_recs_min, n_recs_max in pos_val_l:
#         val_freqs[val] = val_freqs.get(val, 0) + 1.0
#     N = float(sum(val_freqs.values()))
#     val_freqs = val_freqs.items()
#     val_freqs.sort()
#
#     ps = []
#     Nc = N
#     for val, freq in val_freqs:
#         ps.append((val, Nc / N))
#         Nc -= freq
#     return dict(ps)
#
#
# def score_peaks_by_seg(bed_d, segmentation, random_perms=100, flank=15,
#                        analysis_id=None):
#     """ to preserve memory, bed_d will be destroyed as it is being loaded into
#         some temporary structures, which are also destroyed in the end
#     """
#     tot = len(bed_d)
#     prev_prog = 5
#     by_segment = {}
#     while bed_d:
#         (chrome, strand), tmpd = bed_d.popitem()
#         for (start_pos, end_pos), value in tmpd.iteritems():
#             assert (value >= 0.0)
#             hit_int, hit_sense, hit_same, hit_anti = segmentation.get_annotation(
#                 chrome, strand, start_pos)
#             if hit_int is None or hit_sense <> 'same':
#                 continue
#             seg = hit_int.segment
#             seg_id = seg.seg_id
#             by_segment.setdefault((seg_id, chrome, strand), []).append(
#                 (start_pos, value))
#
#         if analysis_id is not None:
#             prog = 25 - (100000 * len(bed_d) / tot) / 5000
#             if prog > prev_prog:
#                 db.analysis_set_status(analysis_id, 'processing %s%%' % prog)
#                 prev_prog = prog
#
#     random.seed(42)
#     ret_d = {}
#     tot = len(by_segment)
#     # process segments always in same order (to obtain same results when sampling randomly)
#     seg_keys = sorted(by_segment.keys())
#     for seg_key in seg_keys:
#         (seg_id, chrome, strand) = seg_key
#         true_hits = None
#         true_hits = by_segment[seg_key]
#         del by_segment[seg_key]
#
#         true_hits.sort()
#         true_hits_extended = extended_sum(true_hits, flank)
#         total_counts = int(sum(val for pos, val in true_hits))
#         seg = segmentation.segments[seg_id]
#
#         # extend true data
#         true_probs = cumulative_probs(true_hits_extended)
#         true_vals = set(true_probs.keys())
#         min_true_val = min(true_vals)
#
#         # randomize and extend
#         all_vals = set(true_vals)
#         rnd_res = []
#         for rnd in xrange(random_perms):
#             rnd_hits = {}
#             for x in range(total_counts):
#                 rnd_pos = seg.random_position()
#                 rnd_hits[rnd_pos] = rnd_hits.get(rnd_pos, 0.0) + 1.0
#             rnd_hits = list(rnd_hits.items())
#             rnd_hits.sort()
#             rnd_hits_extended = extended_sum(rnd_hits, flank)
#             rnd_probs = cumulative_probs(rnd_hits_extended)
#             rnd_res.append(rnd_probs)
#             all_vals.update(set(rnd_probs.keys()))
#         all_vals = [val for val in all_vals if val >= min_true_val]
#         all_vals.sort(reverse=True)
#
#         # by decreasing values (counts)
#         prev_p = [0.0 for i in xrange(random_perms)]
#         p_avg = {}
#         p_std = {}
#         for val in all_vals:
#             rnd_pvals = []
#             for i in xrange(random_perms):
#                 p = rnd_res[i].get(val, prev_p[i])
#                 prev_p[i] = p
#                 rnd_pvals.append(p)
#             if val in true_vals:
#                 p_avg[val] = sum(rnd_pvals) / float(random_perms)
#                 p_std[val] = numpy.std(rnd_pvals)
#
#         # calculate fdr
#         ret_tmpl = ret_d.setdefault((chrome, strand), [])
#         for (start_pos, value), (
#         ext_start_pos, ext_value, n_recs_min, n_recs_max) in zip(true_hits,
#                                                                  true_hits_extended):
#             assert (start_pos == ext_start_pos)
#             p_true = true_probs[ext_value]
#             fdr = (p_avg[ext_value] + p_std[ext_value]) / p_true
#             ret_tmpl.append((start_pos, value, ext_value, n_recs_min,
#                              n_recs_max, p_true, fdr, seg.seg_type,
#                              seg.biotype, seg.gene_name))
#
#         if analysis_id is not None:
#             prog = 60 - (100000 * len(by_segment) / tot) / 2857
#             if prog > prev_prog:
#                 db.analysis_set_status(analysis_id, 'processing %s%%' % prog)
#                 prev_prog = prog
#     # sort all records, on each (chromosome, strand) separately
#     for (chrome, strand), tmpl in ret_d.iteritems():
#         tmpl.sort()
#     return ret_d
#
#
# def join_adjacent(lst):
#     lst.sort()
#     if len(lst) <= 1:
#         return lst
#     new_lst = []
#     _, pe = lst[0]
#     cs, ce = lst[0]
#     for s, e in lst[1:]:
#         if s - pe == 1:
#             ce = e
#         else:
#             new_lst.append((cs, ce))
#             cs, ce = s, e
#         pe = e
#     new_lst.append((cs, ce))
#     return new_lst
#
#
# def score_peaks_as_one(bed_d, segmentation, random_perms=100, flank=15,
#                        analysis_id=None):
#     """ to preserve memory, bed_d will be destroyed as it is being loaded into
#         some temporary structures, which are also destroyed in the end
#     """
#     tot = len(bed_d)
#     prev_prog = 5
#     by_gene = {}
#     while bed_d:
#         (chrome, strand), tmpd = bed_d.popitem()
#         for (start_pos, end_pos), value in tmpd.iteritems():
#             assert (value >= 0.0)
#             hit_int, hit_sense, hit_same, hit_anti = segmentation.get_annotation(
#                 chrome, strand, start_pos)
#             if hit_int is None or hit_sense <> 'same':
#                 continue
#             seg = hit_int.segment
#             if seg.seg_type == 'inter' or seg.seg_type == 'telo':
#                 reg_name = ('non-gene', seg.seg_id)
#             else:
#                 reg_name = ('gene', seg.gene_name)
#             by_gene.setdefault((reg_name, chrome, strand, seg), []).append(
#                 (start_pos, value))
#
#         if analysis_id is not None:
#             prog = 25 - (100000 * len(bed_d) / tot) / 5000
#             if prog > prev_prog:
#                 db.analysis_set_status(analysis_id, 'processing %s%%' % prog)
#                 prev_prog = prog
#
#     random.seed(42)
#     ret_d = {}
#     tot = len(by_gene)
#     progi = 0
#     gene_keys = sorted(by_gene.keys())
#     for gene_key in gene_keys:
#         ((reg_type, reg_name), chrome, strand, orig_seg) = gene_key
#         true_hits = None
#         true_hits = by_gene[gene_key]
#         del by_gene[gene_key]
#
#         true_hits.sort()
#         true_hits_extended = extended_sum(true_hits, flank)
#         total_counts = int(sum(val for pos, val in true_hits))
#
#         # merge all gene segments into one collection of intervals
#         if reg_type == 'non-gene':
#             seg = orig_seg  # segmentation.segments[reg_name]
#         else:
#             seg_chrome, seg_strand, seg_list = segmentation.segments_by_gene[
#                 reg_name]
#             assert (chrome == seg_chrome)
#             assert (strand == seg_strand)
#             seg = annotation.Segment()
#             seg.seg_type = 'ORF'
#             #            print "################################################"
#             #            print "merging", chrome, strand, reg_name
#             for s in seg_list:
#                 #                print "\ta:", s._start_length
#                 #                print "\tb:", s.seg_type, s.intervals
#                 seg.intervals.extend(s.intervals)
#                 seg.intervals_nonconstitutive.extend(
#                     s.intervals_nonconstitutive)
#
#             seg.intervals = join_adjacent(seg.intervals)
#             seg.intervals_nonconstitutive = join_adjacent(
#                 seg.intervals_nonconstitutive)
#             seg._calc_derived_fields()
#         # print "merged"
#         #            print "\ta:", seg._start_length
#         #            print "\tb:", seg.seg_type, seg.intervals
#         #            print
#
#         # extend true data
#         true_probs = cumulative_probs(true_hits_extended)
#         true_vals = set(true_probs.keys())
#         min_true_val = min(true_vals)
#
#         # randomize and extend
#         all_vals = set(true_vals)
#         rnd_res = []
#         for rnd in xrange(random_perms):
#             rnd_hits = {}
#             for x in range(total_counts):
#                 rnd_pos = seg.random_position()
#                 rnd_hits[rnd_pos] = rnd_hits.get(rnd_pos, 0.0) + 1.0
#             rnd_hits = list(rnd_hits.items())
#             rnd_hits.sort()
#             rnd_hits_extended = extended_sum(rnd_hits, flank)
#             rnd_probs = cumulative_probs(rnd_hits_extended)
#             rnd_res.append(rnd_probs)
#             all_vals.update(set(rnd_probs.keys()))
#         all_vals = [val for val in all_vals if val >= min_true_val]
#         all_vals.sort(reverse=True)
#
#         # by decreasing values (counts)
#         prev_p = [0.0 for i in xrange(random_perms)]
#         p_avg = {}
#         p_std = {}
#         for val in all_vals:
#             rnd_pvals = []
#             for i in xrange(random_perms):
#                 p = rnd_res[i].get(val, prev_p[i])
#                 prev_p[i] = p
#                 rnd_pvals.append(p)
#             if val in true_vals:
#                 p_avg[val] = sum(rnd_pvals) / float(random_perms)
#                 p_std[val] = numpy.std(rnd_pvals)
#
#         # calculate fdr
#         ret_tmpl = ret_d.setdefault((chrome, strand), [])
#         for (start_pos, value), (
#         ext_start_pos, ext_value, n_recs_min, n_recs_max) in zip(true_hits,
#                                                                  true_hits_extended):
#             assert (start_pos == ext_start_pos)
#             p_true = true_probs[ext_value]
#             fdr = (p_avg[ext_value] + p_std[ext_value]) / p_true
#             ret_tmpl.append((start_pos, value, ext_value, n_recs_min,
#                              n_recs_max, p_true, fdr, orig_seg.seg_type,
#                              orig_seg.biotype, orig_seg.gene_name))
#
#         if analysis_id is not None:
#             progi += 1
#             prog = 60 - (100000 * progi / tot) / 2857
#             if prog > prev_prog:
#                 db.analysis_set_status(analysis_id, 'processing %s%%' % prog)
#                 prev_prog = prog
#
#     # sort all records, on each (chromosome, strand) separately
#     for (chrome, strand), tmpl in ret_d.iteritems():
#         tmpl.sort()
#     return ret_d
#
#
# def cluster_peaks(peaks, flank):
#     clusters = {}
#     for (chrome, strand), tmpl in peaks.iteritems():
#         if not tmpl:
#             continue
#         tmpd_clusters = {}
#         cl_start, cl_value, cl_recs_min, cl_recs_max = tmpl[0]
#         cl_end = cl_start
#         for start_pos, value, n_recs_min, n_recs_max in tmpl[1:]:
#             if start_pos - cl_end > flank:
#                 # close previous cluster
#                 # !!!assert((cl_start, cl_end+1) not in tmpd_clusters)
#                 tmpd_clusters[(cl_start, cl_end + 1)] = cl_value
#                 # start new cluster
#                 cl_start = start_pos
#                 cl_end = start_pos
#                 cl_value = value
#             else:
#                 # update (end of) current cluster
#                 cl_end = start_pos
#                 cl_value += value
#         tmpd_clusters[(cl_start, cl_end + 1)] = cl_value
#         clusters[(chrome, strand)] = tmpd_clusters
#     return clusters
#
#
# def cluster_peaks_extend(peaks_lowFDR, peaks_all, flank):
#     clusters = {}
#     for (chrome, strand), hits_lowFDR in peaks_lowFDR.iteritems():
#         if not hits_lowFDR:
#             continue
#         tmpl_clusters = []
#         cl_start, _, cl_recs_min, cl_recs_max = hits_lowFDR[0]
#         cl_start = cl_recs_min
#         cl_end = cl_recs_max
#         for _, _, n_recs_min, n_recs_max in hits_lowFDR[1:]:
#             if n_recs_min - cl_end > flank:
#                 # close previous cluster
#                 tmpl_clusters.append((cl_start, cl_end))
#                 # start new cluster
#                 cl_start = n_recs_min
#                 cl_end = n_recs_max
#             else:
#                 # update (end of) current cluster
#                 cl_end = n_recs_max
#         tmpl_clusters.append((cl_start, cl_end))
#
#         # iterate over list of clusters and peaks, and sum counts within each cluster
#         clusters_with_counts = {}
#         hits_all = peaks_all.get((chrome, strand), [])
#         L_peaks = len(hits_all)
#         cur_peak = 0
#         for cl_start, cl_end in tmpl_clusters:
#             while cur_peak < L_peaks and hits_all[cur_peak][0] < cl_start:
#                 cur_peak += 1
#
#             if cur_peak >= L_peaks:
#                 clusters_with_counts[(cl_start, cl_end + 1)] = 0.0
#                 break
#
#             cur_sum = 0.0
#             while cur_peak < L_peaks and hits_all[cur_peak][0] <= cl_end:
#                 cur_sum += hits_all[cur_peak][1]
#                 cur_peak += 1
#             clusters_with_counts[(cl_start, cl_end + 1)] = cur_sum
#         clusters[(chrome, strand)] = clusters_with_counts
#     return clusters
#
#
# class CountStats:
#     def __init__(self):
#         self.same_count = 0.0
#         self.anti_count = 0.0
#         self.same_sites = 0.0
#         self.anti_sites = 0.0
#
#     def add_same(self, count, sites=1):
#         self.same_count += count
#         self.same_sites += sites
#
#     def add_anti(self, count, sites=1):
#         self.anti_count += count
#         self.anti_sites += sites
#
#
# class Gen_base:
#     def __init__(self, out_filename, custom_annot,
#                  write_custom_annot_line=False):
#         self.out_filename = out_filename
#         self.f = iCount.nc_open(self.out_filename, "wt")
#         self.description = self.name = os.path.basename(out_filename)
#         self.records = []
#         Lcustom_annot = {}
#         Lcustom_annot.update(custom_annot)
#         if write_custom_annot_line:
#             self.f.write("%s\n" % ", ".join(
#                 ["%s=%s" % (a, v) for a, v in Lcustom_annot.iteritems()]))
#
#     def close(self):
#         self.write_records()
#         self.f.close()
#
#
# class Gen_sumsegt(Gen_base):
#     header = [
#         "chromosome", "strand", "start_position", "end_position", "segment_id",
#         "segment_gene_name", "segment_type", "segment_biotype",
#         "segment_length",
#         "same_count_sum", "same_count_density",
#         "anti_count_sum", "anti_count_density",
#         "same_sites", "same_sites_density",
#         "anti_sites", "anti_sites_density"
#     ]
#
#     def __init__(self, out_filename, custom_annot, segmentation):
#         Gen_base.__init__(self, out_filename, custom_annot,
#                           write_custom_annot_line=True)
#         self.f.write("%s\n" % "\t".join(self.header))
#         self.stats_by_segment = {}
#         self.segmentation = segmentation
#
#     def update_from_counts(self, chrome, hit_pos, strand, count, pann):
#         k = (pann.seg_id, chrome)
#         if k in self.stats_by_segment:
#             cs = self.stats_by_segment[k]
#         else:
#             cs = CountStats()
#             self.stats_by_segment[k] = cs
#         if pann.seg_strand == "same":
#             cs.add_same(count)
#         else:
#             cs.add_anti(count)
#
#     def write_records(self):
#         cs_by_seg_id = {}
#         # write records on undefined genome segments at beginning of file
#         # there should be no such segments, except for those chromosomes that are not in biomart
#         for (seg_id, chrome), cs in self.stats_by_segment.iteritems():
#             if seg_id == "none":
#                 r = (chrome, '?', '?', '?', seg_id, 'undefined', 'undefined',
#                      'undefined', '?',
#                      "%g" % cs.same_count, "?",
#                      "%g" % cs.anti_count, "?",
#                      "%g" % cs.same_sites, "?",
#                      "%g" % cs.anti_sites, "?",
#                      )
#                 assert (len(r) == len(self.header))
#                 self.f.write("%s\n" % "\t".join([str(x) for x in r]))
#             else:
#                 # !!!assert(seg_id not in cs_by_seg_id)
#                 cs_by_seg_id[seg_id] = (chrome, cs)
#
#         # write records on "valid" segments
#         for seg_id in self.segmentation.segments_order:
#             sr = self.segmentation.segments[seg_id]
#             if sr.strand == '+':
#                 start_position = sr.min_pos
#                 end_position = sr.max_pos
#             else:
#                 start_position = sr.max_pos
#                 end_position = sr.min_pos
#
#             if seg_id in cs_by_seg_id:
#                 chrome, cs = cs_by_seg_id[seg_id]
#                 assert (sr.chromosome == chrome and sr.seg_id == seg_id)
#                 seg_len = sr.total_len
#
#                 # same sense statistics
#                 same_count_density = cs.same_count / float(seg_len)
#                 same_sites_density = cs.same_sites / float(seg_len)
#
#                 # anti sense statistics
#                 anti_count_density = cs.anti_count / float(seg_len)
#                 anti_sites_density = cs.anti_sites / float(seg_len)
#
#                 r = (
#                     sr.chromosome, sr.strand, str(start_position),
#                     str(end_position), seg_id, sr.gene_name, sr.seg_type,
#                     sr.biotype, seg_len,
#                     "%g" % cs.same_count, "%g" % same_count_density,
#                     "%g" % cs.anti_count, "%g" % anti_count_density,
#                     "%g" % cs.same_sites, "%g" % same_sites_density,
#                     "%g" % cs.anti_sites, "%g" % anti_sites_density,
#                 )
#             else:
#                 r = (sr.chromosome, sr.strand, str(start_position),
#                      str(end_position), seg_id, sr.gene_name, sr.seg_type,
#                      sr.biotype, sr.total_len) + ('',) * (len(self.header) - 9)
#
#             assert (len(r) == len(self.header))
#             self.f.write("%s\n" % "\t".join([str(x) for x in r]))
#
#
# class Gen_sumgene(Gen_base):
#     types_for_sumgene = ["total", "ncRNA", "ORF", "5UTR", "3UTR", "intron"]
#     header = ["chromosome", "strand", "start_position", "end_position",
#               "gene_name", "gene_segments_types", "gene_segments_biotypes"]
#
#     for seg_type in types_for_sumgene:
#         header.extend(
#             ["%s_%s" % (seg_type, x) for x in [
#                 "same_count_sum", "same_count_density",
#                 "anti_count_sum", "anti_count_density",
#                 "same_sites", "same_sites_density",
#                 "anti_sites", "anti_sites_density",
#             ]
#              ]
#         )
#
#     def __init__(self, out_filename, custom_annot, segmentation):
#         Gen_base.__init__(self, out_filename, custom_annot,
#                           write_custom_annot_line=True)
#         self.f.write("%s\n" % "\t".join(self.header))
#         self.stats_by_gene = {}
#         self.segmentation = segmentation
#
#     def update_from_counts(self, chrome, hit_pos, strand, count, pann):
#         if pann.seg_type not in self.types_for_sumgene:
#             return
#
#         k = (pann.seg_gene_name, chrome)
#         if k in self.stats_by_gene:
#             tmpd = self.stats_by_gene[k]
#         else:
#             self.stats_by_gene[k] = tmpd = {}
#
#         if pann.seg_type in tmpd:
#             cs, _ = tmpd[pann.seg_type]
#         else:
#             cs = CountStats()
#             tmpd[pann.seg_type] = cs, pann.seg
#         if pann.seg_strand == "same":
#             cs.add_same(count)
#         else:
#             cs.add_anti(count)
#
#     def write_records(self):
#         r_by_gene = {}
#         # write records on undefined genome segments at beginning of file
#         # there should be no such segments, except for those chromosomes that are not in biomart
#         for (gene_name, chrome), tmpd in self.stats_by_gene.iteritems():
#             if gene_name == 'undefined':
#                 rr = []
#                 cs_sum = CountStats()
#                 for seg_type in self.types_for_sumgene:
#                     if seg_type in tmpd:
#                         cs, pann_seg = tmpd[seg_type]
#                         rr.extend([
#                             "%g" % cs.same_count, "?",
#                             "%g" % cs.anti_count, "?",
#                             "%g" % cs.same_sites, "?",
#                             "%g" % cs.anti_sites, "?",
#                         ])
#                         cs_sum.add_same(cs.same_count, cs.same_sites)
#                         cs_sum.add_anti(cs.anti_count, cs.anti_sites)
#                     else:
#                         rr.extend(['0'] * 8)
#                 r = [chrome, '?', '?', '?', 'undefined', 'undefined',
#                      'undefined']
#                 r.extend([
#                     "%g" % cs_sum.same_count, "?",
#                     "%g" % cs_sum.anti_count, "?",
#                     "%g" % cs_sum.same_sites, "?",
#                     "%g" % cs_sum.anti_sites, "?",
#                 ])
#                 r.extend(rr)
#                 assert (len(r) == len(self.header))
#                 self.f.write("%s\n" % "\t".join([str(x) for x in r]))
#             else:
#                 # !!!assert(gene_name not in r_by_gene)
#                 rr = []
#                 cs_sum = CountStats()
#                 cs_den = CountStats()
#                 assert (self.types_for_sumgene[0] == "total")
#                 for seg_type in self.types_for_sumgene[
#                                 1:]:  # skip total, since it is calculated last
#                     if seg_type in tmpd:
#                         cs, pann_seg = tmpd[seg_type]
#                         seg_len = pann_seg.total_len
#
#                         # same sense statistics
#                         same_count_density = cs.same_count / float(seg_len)
#                         same_sites_density = cs.same_sites / float(seg_len)
#
#                         # anti sense statistics
#                         anti_count_density = cs.anti_count / float(seg_len)
#                         anti_sites_density = cs.anti_sites / float(seg_len)
#
#                         rr.extend([
#                             "%g" % cs.same_count, "%g" % same_count_density,
#                             "%g" % cs.anti_count, "%g" % anti_count_density,
#                             "%g" % cs.same_sites, "%g" % same_sites_density,
#                             "%g" % cs.anti_sites, "%g" % anti_sites_density,
#                         ])
#
#                         # total sum
#                         cs_sum.add_same(cs.same_count, cs.same_sites)
#                         cs_sum.add_anti(cs.anti_count, cs.anti_sites)
#                         # total density
#                         cs_den.add_same(same_count_density, same_sites_density)
#                         cs_den.add_anti(anti_count_density, anti_sites_density)
#                     else:
#                         rr.extend(['0'] * 8)
#                 r = [chrome, '?', '?', '?', gene_name, '?', '?']
#                 r.extend([
#                     "%g" % cs_sum.same_count, "%g" % cs_den.same_count,
#                     "%g" % cs_sum.anti_count, "%g" % cs_den.anti_count,
#                     "%g" % cs_sum.same_sites, "%g" % cs_den.same_sites,
#                     "%g" % cs_sum.anti_sites, "%g" % cs_den.anti_sites,
#                 ])
#                 r.extend(rr)
#                 # !!!assert(gene_name not in r_by_gene)
#                 r_by_gene[gene_name] = r
#
#         for gene_name in self.segmentation.genes_order:
#             chrome, strand, segments = self.segmentation.segments_by_gene[
#                 gene_name]
#             if strand == '+':
#                 start_position = min([sr.min_pos for sr in segments])
#                 end_position = max([sr.max_pos for sr in segments])
#             else:
#                 start_position = max([sr.max_pos for sr in segments])
#                 end_position = min([sr.min_pos for sr in segments])
#
#             seg_types = sorted(list(set(sr.seg_type for sr in segments)))
#             seg_types = ",".join(seg_types) if seg_types else "undefined"
#             seg_biotypes = sorted(list(set(sr.biotype for sr in segments)))
#             seg_biotypes = ",".join(
#                 seg_biotypes) if seg_biotypes else "undefined"
#             r = r_by_gene.get(gene_name, [chrome, strand, str(start_position),
#                                           str(end_position), gene_name,
#                                           seg_types, seg_biotypes] + [''] * (
#                               len(self.header) - 7))
#             r[1] = strand
#             r[2] = str(start_position)
#             r[3] = str(end_position)
#             r[5] = seg_types
#             r[6] = seg_biotypes
#             assert (r[0] == chrome and r[1] == strand and r[4] == gene_name)
#             assert (len(r) == len(self.header))
#             self.f.write("%s\n" % "\t".join(r))
#
#
# class Gen_sumtype(Gen_base):
#     header = [
#         "segment_type",
#         "same_count_sum", "%", "anti_count_sum", "%",
#         "same_sites", "%", "anti_sites", "%"
#     ]
#
#     header_density_enrichment = [
#         "segment_type", "segment_length", "%",
#         "same_count_sum_density_enrichment",
#         "anti_count_sum_density_enrichment",
#         "same_sites_density_enrichment", "anti_sites_density_enrichment"
#     ]
#
#     def __init__(self, out_filename, custom_annot, segmentation,
#                  include_biotype=False):
#         Gen_base.__init__(self, out_filename, custom_annot,
#                           write_custom_annot_line=True)
#         self.f.write("%s\n" % "\t".join(self.header))
#         self.stats_by_type = {}
#         self.include_biotype = include_biotype
#         self.segmentation = segmentation
#
#     def update_from_counts(self, chrome, hit_pos, strand, count, pann):
#         if self.include_biotype:
#             if pann.seg_biotype and pann.seg_type <> 'undefined':
#                 k = "%s, %s" % (pann.seg_type, pann.seg_biotype)
#             else:
#                 k = "%s" % (pann.seg_type)
#         else:
#             k = pann.seg_type
#         if k in self.stats_by_type:
#             cs = self.stats_by_type[k]
#         else:
#             cs = CountStats()
#             self.stats_by_type[k] = cs
#         if pann.seg_strand == "same":
#             cs.add_same(count)
#         else:
#             cs.add_anti(count)
#
#     def write_records(self):
#         self.records = []
#         self.records_density_enrichment = []
#         same_count = 0.0
#         same_sites = 0.0
#         anti_count = 0.0
#         anti_sites = 0.0
#         for seg_type, cs in self.stats_by_type.iteritems():
#             same_count += cs.same_count
#             same_sites += cs.same_sites
#             anti_count += cs.anti_count
#             anti_sites += cs.anti_sites
#         tot_count = same_count + anti_count
#         tot_sites = same_sites + anti_sites
#         if tot_count == 0.0: tot_count = 1.0
#         if tot_sites == 0.0: tot_sites = 1.0
#         if self.include_biotype:
#             type_keys = self.segmentation.total_len_biotype.keys()
#             seg_total_length = sum(
#                 self.segmentation.total_len_biotype.values())
#         else:
#             type_keys = self.segmentation.total_len_seg_type.keys()
#             seg_total_length = sum(
#                 self.segmentation.total_len_seg_type.values())
#         seg_total_per_len_sum = 0.0
#         for seg_type in type_keys:  # self.stats_by_type.iterkeys():
#             cs = self.stats_by_type.get(seg_type, CountStats())
#             if self.include_biotype:
#                 seg_length = self.segmentation.total_len_biotype.get(seg_type,
#                                                                      0)
#             else:
#                 seg_length = self.segmentation.total_len_seg_type.get(seg_type,
#                                                                       0)
#             per_len = 100.0 * seg_length / seg_total_length
#             seg_total_per_len_sum += per_len
#
#             per_same_count = 100.0 * cs.same_count / tot_count
#             per_anti_count = 100.0 * cs.anti_count / tot_count
#             per_same_sites = 100.0 * cs.same_sites / tot_sites
#             per_anti_sites = 100.0 * cs.anti_sites / tot_sites
#
#             self.records.append((seg_type,
#                                  "%g" % cs.same_count,
#                                  "%0.2f%%" % per_same_count,
#                                  "%g" % cs.anti_count,
#                                  "%0.2f%%" % per_anti_count,
#                                  "%g" % cs.same_sites,
#                                  "%0.2f%%" % per_same_sites,
#                                  "%g" % cs.anti_sites,
#                                  "%0.2f%%" % per_anti_sites,
#                                  ))
#
#             if per_len > 0.0:
#                 self.records_density_enrichment.append((seg_type,
#                                                         "%s" % seg_length,
#                                                         "%g%%" % per_len,
#                                                         "%g" % (
#                                                         per_same_count / per_len),
#                                                         "%g" % (
#                                                         per_anti_count / per_len),
#                                                         "%g" % (
#                                                         per_same_sites / per_len),
#                                                         "%g" % (
#                                                         per_anti_sites / per_len),
#                                                         ))
#             else:
#                 self.records_density_enrichment.append(
#                     (seg_type, "%s" % seg_length, "%g%%" % per_len,) + (
#                     '?',) * 4)
#
#         self.records.sort()
#         self.records.append((
#             "SUM",
#             "%g" % same_count, "%0.2f%%" % (100.0 * same_count / tot_count),
#             "%g" % anti_count, "%0.2f%%" % (100.0 * anti_count / tot_count),
#             "%g" % same_sites, "%0.2f%%" % (100.0 * same_sites / tot_sites),
#             "%g" % anti_sites, "%0.2f%%" % (100.0 * anti_sites / tot_sites),
#         ))
#         self.records_density_enrichment.sort()
#         self.records_density_enrichment.append((
#             "SUM", "%s" % seg_total_length, "%g%%" % (seg_total_per_len_sum),
#         ))
#         for r in self.records:
#             self.f.write("%s\n" % "\t".join(r))
#         self.f.write("\n")
#         self.f.write("Density enrichments:\n")
#         self.f.write("%s\n" % "\t".join(self.header_density_enrichment))
#         for r in self.records_density_enrichment:
#             self.f.write("%s\n" % "\t".join(r))
#
#
# def process_file(in_fname, mapped_to, annot_ver, out_fname_pref, random_perms,
#                  flank, fdr_th, regions, segmentation=None, analysis_id=None):
#     print
#     "*********************************"
#     print
#     time.ctime()
#     print
#     "searching peaks in file:", in_fname
#     print
#     "saving to:", out_fname_pref
#     print
#     "parameters:"
#     print
#     "\trandom permutations:", random_perms
#     print
#     "\tflanking:", flank
#     print
#     "\tfdr threshold:", fdr_th
#     print
#     "\tregions:", regions
#     sys.stdout.flush()
#
#     if analysis_id is not None:
#         db.analysis_set_status(analysis_id, 'processing started')
#         prev_prog = 0
#
#     bed_d, bed_h = BED.load_bedGraph(in_fname, return_tracks_descriptions=True)
#     assert (len(bed_h) == 1)
#     bed_h = bed_h[0]
#     try:
#         bed_res_type = bed_h.split('res_type=')[1].split(' ')[0].strip('"')
#     except:
#         bed_res_type = 'unknown_counts'
#     bed_res_type_orig = bed_res_type
#     bed_res_type_desc = beds.desc_by_result_type.get(bed_res_type, None)
#     if bed_res_type_desc is not None:
#         bed_res_type = bed_res_type_desc[-1]
#     print
#     "input res_type: %s" % bed_res_type
#     sys.stdout.flush()
#
#     if segmentation is None or segmentation.assembly <> mapped_to or segmentation.annotation_version <> annot_ver:
#         segmentation = iCount.genomes.annotation.Segmentation(
#             assembly=mapped_to, annotation_version=annot_ver)
#         print
#         "segmentation loaded"
#         sys.stdout.flush()
#
#     if analysis_id is not None:
#         db.analysis_set_status(analysis_id, 'processing 5%')
#         prev_prog = 5
#
#     if regions == 'separately':
#         peaks = score_peaks_by_seg(bed_d, segmentation, random_perms, flank,
#                                    analysis_id=analysis_id)
#     else:
#         peaks = score_peaks_as_one(bed_d, segmentation, random_perms, flank,
#                                    analysis_id=analysis_id)
#
#     all_chromes = list(set([chrome for chrome, strand in peaks.keys()]))
#     all_chromes.sort()
#     print
#     "peaks identified"
#     sys.stdout.flush()
#
#     custom_annot = {
#         "mapped_to": mapped_to, "annot_ver": annot_ver,
#         "random_perms": random_perms, "flank": flank,
#         "fdr_th": fdr_th,
#         "regions": regions,
#         "input_file": os.path.basename(in_fname),
#         "res_type": bed_res_type_orig,
#         "res_type_values": bed_res_type
#     }
#
#     # summary files
#     out_sumgene_gen = Gen_sumgene(out_fname_pref + "_sumgene.tab.gz",
#                                   custom_annot, segmentation)
#     out_sumsegt_gen = Gen_sumsegt(out_fname_pref + "_sumsegt.tab.gz",
#                                   custom_annot, segmentation)
#     out_sumtype_gen = Gen_sumtype(out_fname_pref + "_sumtype.tab.gz",
#                                   custom_annot, segmentation,
#                                   include_biotype=False)
#     out_sumbtype_gen = Gen_sumtype(out_fname_pref + "_sumbtype.tab.gz",
#                                    custom_annot, segmentation,
#                                    include_biotype=True)
#
#     # annotated positions
#     # generate lowFDR.bed records - report only positions with FDR < 0.05
#     # _scores.tab
#     fout = iCount.nc_open(out_fname_pref + "_scores.tab.gz", "wt")
#     fout.write("%s\n" % "\t".join(
#         ['chromosome', 'hit_position', 'strand', 'seg_type', 'seg_biotype',
#          'gene_name', 'value [%s]' % bed_res_type,
#          'extended_value_%snt [%s]' % (flank, bed_res_type), 'p.value',
#          'FDR']))
#     peaks_lowFDR = {}
#     bed_lowFDR = []
#     for i, chrome in enumerate(all_chromes):
#         recs_out = []
#         recs_out_bed = []
#         for strand in ['+', '-']:
#             recs_peaks = peaks.get((chrome, strand), [])
#             tmpl_peaks_low = peaks_lowFDR.setdefault((chrome, strand), [])
#             for start_pos, value, ext_value, n_recs_min, n_recs_max, p_true, fdr, seg_type, seg_biotype, seg_gene_name in recs_peaks:
#                 recs_out.append((chrome, start_pos, strand, seg_type,
#                                  seg_biotype, seg_gene_name, value, ext_value,
#                                  p_true, fdr))
#                 if fdr < fdr_th:
#                     # annotate position and add all summary statistics
#                     pos_annot = beds.AnnotatePosition(segmentation, chrome,
#                                                       strand, start_pos, None)
#                     out_sumgene_gen.update_from_counts(chrome, start_pos,
#                                                        strand, value,
#                                                        pos_annot)
#                     out_sumsegt_gen.update_from_counts(chrome, start_pos,
#                                                        strand, value,
#                                                        pos_annot)
#                     out_sumtype_gen.update_from_counts(chrome, start_pos,
#                                                        strand, value,
#                                                        pos_annot)
#                     out_sumbtype_gen.update_from_counts(chrome, start_pos,
#                                                         strand, value,
#                                                         pos_annot)
#
#                     # add for output in bedGraph file
#                     tmpl_peaks_low.append(
#                         (start_pos, value, n_recs_min, n_recs_max))
#                     if strand == '-':
#                         value = -value
#                     recs_out_bed.append(
#                         (chrome, start_pos, start_pos + 1, value))
#
#         recs_out.sort()
#         for chrome, start_pos, strand, seg_type, seg_biotype, seg_gene_name, value, ext_value, p_true, fdr in recs_out:
#             olst = [chrome, start_pos, strand, seg_type, seg_biotype,
#                     seg_gene_name, "%s%g" % (strand, value),
#                     "%s%g" % (strand, ext_value), "%g" % (p_true),
#                     "%g" % (fdr)]
#             fout.write("%s\n" % "\t".join([str(x) for x in olst]))
#         recs_out_bed.sort()
#         bed_lowFDR.extend(recs_out_bed)
#         del recs_out
#         del recs_out_bed
#
#         if analysis_id is not None:
#             prog = 60 + (100000 * (i + 1) / len(all_chromes)) / 4000
#             if prog > prev_prog:
#                 db.analysis_set_status(analysis_id, 'processing %s%%' % prog)
#                 prev_prog = prog
#
#     fout.close()
#     print
#     "peaks scores stored to:", out_fname_pref + "_scores.tab.gz"
#     sys.stdout.flush()
#     # save summaries
#     out_sumgene_gen.close()
#     out_sumsegt_gen.close()
#     out_sumtype_gen.close()
#     out_sumbtype_gen.close()
#
#     # _lowFDR.bed
#     name = os.path.basename(out_fname_pref)
#     desc = "peaks in %s, with %s permutations, %s nt neighborhood, FDR < %g, regions %s" % (
#     os.path.basename(in_fname), random_perms, flank, fdr_th, regions)
#     BED.save_bedGraph(out_fname_pref + "_lowFDR.bed.gz", bed_lowFDR,
#                       name + "_lowFDR.bed.gz", desc, mapped_to)
#     del bed_lowFDR
#     if analysis_id is not None:
#         db.analysis_set_status(analysis_id, 'processing 90%')
#     print
#     "peaks stored to:", out_fname_pref + "_lowFDR.bed.gz"
#     sys.stdout.flush()
#
#     # clustered - cluster neighboring peaks, within +-15nt, report value (cDNA, tc,...) sum for region covered by cluster
#     clusters = cluster_peaks(peaks_lowFDR, flank)
#     desc = "peaks in %s, with %s permutations, FDR < %g, regions %s, clustered within %s nt neighborhood" % (
#     os.path.basename(in_fname), random_perms, fdr_th, regions, flank)
#     BED.save_bedGraph(out_fname_pref + "_lowFDR_clusters.bed.gz", clusters,
#                       name + "_lowFDR_clusters.bed.gz", desc, mapped_to)
#     del clusters
#     if analysis_id is not None:
#         db.analysis_set_status(analysis_id, 'processing 95%')
#     print
#     "peaks clustered and stored to:", out_fname_pref + "_lowFDR_clusters.bed.gz"
#     sys.stdout.flush()
#
#     # extended - extend clusters with neighboring non-lowFDR peaks, join clusters if 15nt or less apart
#     clusters_extended = cluster_peaks_extend(peaks_lowFDR, peaks, flank)
#     desc = "peaks in %s, with %s permutations, FDR < %g, regions %s, clustered within %s nt neighborhood, extended to neighboring high FDR peaks or clusters %s nt or less apart" % (
#     os.path.basename(in_fname), random_perms, fdr_th, regions, flank, flank)
#     BED.save_bedGraph(out_fname_pref + "_lowFDR_clusters_ext.bed.gz",
#                       clusters_extended, name + "_lowFDR_clusters_ext.bed.gz",
#                       desc, mapped_to)
#     del clusters_extended
#     print
#     "extended clustered peaks stored to:", out_fname_pref + "_lowFDR_clusters_ext.bed.gz"
#     sys.stdout.flush()
#     if analysis_id is not None:
#         files = [
#             ('tab-peak_scores',
#              iCount.strip_storage_prefix(out_fname_pref + "_scores.tab.gz"),
#              db.filesize_nicefmt(out_fname_pref + "_scores.tab.gz")
#              ),
#             ('tab-sumgene',
#              iCount.strip_storage_prefix(out_fname_pref + "_sumgene.tab.gz"),
#              db.filesize_nicefmt(out_fname_pref + "_sumgene.tab.gz")
#              ),
#             ('tab-sumsegt',
#              iCount.strip_storage_prefix(out_fname_pref + "_sumsegt.tab.gz"),
#              db.filesize_nicefmt(out_fname_pref + "_sumsegt.tab.gz")
#              ),
#             ('tab-sumtype',
#              iCount.strip_storage_prefix(out_fname_pref + "_sumtype.tab.gz"),
#              db.filesize_nicefmt(out_fname_pref + "_sumtype.tab.gz")
#              ),
#             ('tab-sumbtype',
#              iCount.strip_storage_prefix(out_fname_pref + "_sumbtype.tab.gz"),
#              db.filesize_nicefmt(out_fname_pref + "_sumbtype.tab.gz")
#              ),
#             ('bedGraph-lowFDR',
#              iCount.strip_storage_prefix(out_fname_pref + "_lowFDR.bed.gz"),
#              db.filesize_nicefmt(out_fname_pref + "_lowFDR.bed.gz")
#              ),
#             ('bed-lowFDR_clusters',
#              iCount.strip_storage_prefix(
#                  out_fname_pref + "_lowFDR_clusters.bed.gz"),
#              db.filesize_nicefmt(out_fname_pref + "_lowFDR_clusters.bed.gz")
#              ),
#             ('bed-lowFDR_clusters_extended',
#              iCount.strip_storage_prefix(
#                  out_fname_pref + "_lowFDR_clusters_ext.bed.gz"),
#              db.filesize_nicefmt(
#                  out_fname_pref + "_lowFDR_clusters_ext.bed.gz")
#              ),
#         ]
#         db.analysis_set_output(analysis_id, files)
#         db.analysis_set_status(analysis_id, 'done')
#         db.analysis_set_last_update(analysis_id)
#     print
#     "DONE"
#     print
#     time.ctime()
#     print
#     "*********************************"
#     print
#     sys.stdout.flush()
#
