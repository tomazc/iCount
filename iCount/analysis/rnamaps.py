# pylint: skip-file
""".. Line to protect from pydocstyle D205, D400.

RNA maps
--------

Distribution of cross-links relative to genomic landmarks.
Read landmark annotation.
Read bedGraph with cross-linked sites.
Count frequency of cross-links relative to annotation features.
Draw histograms and save tables with distributions.

"""

import logging

import iCount

LOGGER = logging.getLogger(__name__)


# description and parameters needed for the analysis
analysis_name = 'rnamaps'
analysis_description_short = 'RNA maps'
analysis_description = 'Distribution of cross-linked sites relative to ' \
                       'genomic regions.'

params_opt = [
    (
        'maxUp', 'int_range', (-400, -5000, -100), False,
        'Maximum distance to upstream neighboring region.'
    ),
    (
        'maxDown', 'int_range', (400, 100, 5000), False,
        'Maximum distance to downstream neighboring region.'
    ),
    (
        'smooth', 'int_range', (5, 0, 30), False,
        'Half-width of Gaussian window used to smooth RNAmap histogram ('
        'before normalization)'
    ),
    (
        'chromosomes', 'str_list', [], False,
        'Consider only cross-links on listed chromosomes.'

    ),
]

params_pos = [
    (
        'annotation', 'GTF', 'in', 1,
        '(input) GTF file with genomic regions (intervals) and landmarks ('
        'individual sites).'
    ),
    (
        'sites', 'bedGraph', 'in', 1,
        '(input) bedGraph with cross-linked sites.'
    ),
    (
        'table', 'tsv', 'out', 1,
        '(output) Tab-delimited file with distributions.'
    ),
    (
        'maps', 'pdf', 'out', 1,
        '(output) Rendering of the distributions (PDF file).'
    ),
]

# import matplotlib
#
# matplotlib.use("Agg")
# import matplotlib.pyplot as plt
#
# matplotlib.rcParams['axes.labelsize'] = 7  # 'small'
# matplotlib.rcParams['axes.titlesize'] = 7
# matplotlib.rcParams['xtick.labelsize'] = 6  # 'small'
# matplotlib.rcParams['ytick.labelsize'] = 6  # 'small'
# matplotlib.rcParams['legend.fontsize'] = 7

# import time
# import sys
# import os
# import stat
# import glob
# import bisect
# import math
#
# import iCount
# from iCount.files import BED
# from iCount import db
# from iCount import beds

# max_flank = 5000
# # region goes from -5000..+5000
# # position zero (0) is at the beginning of the second type (for given first-second region type ma
# p)
#
# # needed to generate RNAmaps and also to generate definition of RNAmaps
# maps_to_gen = [
#     # map_name, left_side_types, right_side_types, inside gene, inside gene, within same gene, mi
# nU, minD
#     # inside gene = indicates whether part of junction inside gene or intergenic
#     # within same gene
#     ('exon-intron', ['ORF', '3UTR', '5UTR'], ['intron'], True, True, True, -30,
#      100),
#     (
#     'intron-exon', ['intron'], ['ORF', '3UTR', '5UTR'], True, True, True, -100,
#     30),
#     ('5-ncRNA', ['intron', 'inter', 'telo'], ['ncRNA'], False, True, False,
#      -100, 20),
#     (
#     'ncRNA-3', ['ncRNA'], ['intron', 'inter', 'telo'], True, False, False, -20,
#     100),
#     ('ORF-3UTR', ['ORF'], ['3UTR'], True, True, True, -50, 100),
#     ('5UTR-ORF', ['5UTR'], ['ORF'], True, True, True, -50, 50),
#     ("inter-5UTR", ['inter', 'telo'], ['5UTR'], False, True, False, -100, 50),
#     ('3UTR-inter', ['3UTR'], ['inter', 'telo'], True, False, False, -100, 100),
#     ('inter-5', ['inter', 'telo'], ['ORF', 'intron', '3UTR', '5UTR', 'ncRNA'],
#      False, True, False, -100, 50),
#     ('3-inter', ['ORF', 'intron', '3UTR', '5UTR', 'ncRNA'], ['inter', 'telo'],
#      True, False, False, -100, 100),
# ]
#
#
# def _fname_RNAmap_folder(aname, aver):
#     return os.path.join(iCount.genomes.storage_path,
#                         aname, "annotation", aver, "RNAmaps")
#
#
# def _fname_RNAmap(aname, aver, map_name):
#     return os.path.join(_fname_RNAmap_folder(aname, aver),
#                         "%s_%s_RNAmap_%s.tab" % (aname, aver, map_name))
#
#
# class RNAmap:
#     def __init__(self, rnamap_fname='', maxU_within=-300, maxD_within=300,
#                  maxU_flanking=-1000, maxD_flanking=1000):
#         # user-defined
#         self.rnamap_fname = rnamap_fname
#         self.maxU_within = maxU_within
#         self.maxD_within = maxD_within
#         self.maxU_flanking = maxU_flanking
#         self.maxD_flanking = maxD_flanking
#
#         # loaded from definition file
#         self.name = ''
#         self.kind = None
#         self.maxU = None
#         self.maxD = None
#         self.descU = ''  # description of upstream region
#         self.descD = ''  # description of downstream region
#         self.junctions = {}  # key: (chrome, strand), item: [ordered list of junction records]
#         self.distrib = {}  # key: position relative to junction, item: number of junctions spanni
# ng position
#         self.jun_pos2gene_order = {}  # key: (chrome, strand): item: {key:pos: item: (order from
# start, order from end)}
#
#         # temporary structures needed to speed-up search
#         self._bisect_juncs_sta = {}  # key: (chrome, strand), item: [ordered list of junction sta
# rts] needed for bisection
#         self._cur_strand = None
#         self._cur_strand_fact = 0
#         self._cur_juncs_same = None
#         self._cur_juncs_anti = None
#         self._cur_bisect_juncs_sta_same = None
#         self._cur_bisect_juncs_sta_anti = None
#
#         if self.rnamap_fname:
#             self.load(self.rnamap_fname)
#
#     def load(self, rnamap_fname):
#         self.rnamap_fname = rnamap_fname
#
#         ## RNAmap definition
#         f = open(self.rnamap_fname, "rt")
#         name, kind, junctions_cn = [x.strip() for x in
#                                     f.readline().rstrip('\n\r').split('\t')]
#         max_u_span, max_d_span, u_span, d_span = [x.strip() for x in
#                                                   f.readline().rstrip(
#                                                       '\n\r').split('\t')]
#         self.descU = f.readline().rstrip('\n\r')
#         self.descD = f.readline().rstrip('\n\r')
#         max_u_span = int(max_u_span)
#         max_d_span = int(max_d_span)
#         junctions_cn = int(junctions_cn)
#         self.name = name
#         self.kind = kind
#         assert u_span in ['within', 'flanking']
#         assert d_span in ['within', 'flanking']
#         if u_span == 'within':
#             self.maxU = self.maxU_within
#         else:
#             self.maxU = self.maxU_flanking
#         self.maxU = max(self.maxU, max_u_span)
#
#         if d_span == 'within':
#             self.maxD = self.maxD_within
#         else:
#             self.maxD = self.maxD_flanking
#         self.maxD = min(self.maxD, max_d_span)
#
#         ## junction coordinates
#         self.junctions = {}  # (chrome, strand)
#         _tmpd_cache = {}
#         _tmp_jp2go = {}  # make a list of junction positions within each gene
#         for cn in range(junctions_cn):
#             chrome, strand, jun_pos, spanU, spanD, associated_gene_name, associated_gene_ID = f.r
# eadline().rstrip(
#                 '\n\r').split('\t')
#             chrome = _tmpd_cache.setdefault(chrome, chrome)
#             strand = _tmpd_cache.setdefault(strand, strand)
#             jun_pos = int(jun_pos)
#             spanU = max(int(spanU), self.maxU)
#             spanD = min(int(spanD), self.maxD)
#             # calculate jun_start and jun_stop on positive strand
#             if strand == '+':
#                 jun_start = jun_pos + spanU
#                 jun_stop = jun_pos + spanD
#             else:
#                 jun_start = jun_pos - spanD
#                 jun_stop = jun_pos - spanU
#             self.junctions.setdefault((chrome, strand), []).append((jun_start,
#                                                                     jun_stop,
#                                                                     chrome,
#                                                                     strand,
#                                                                     jun_pos,
#                                                                     spanU,
#                                                                     spanD,
#                                                                     associated_gene_name,
#                                                                     associated_gene_ID))
#             _tmp_jp2go.setdefault((chrome, strand), {}).setdefault(
#                 associated_gene_ID, []).append(jun_pos)
#
#         self.jun_pos2gene_order = {}
#         for (chrome, strand), by_gene_ID in _tmp_jp2go.iteritems():
#             tmpd = {}
#             for associated_gene_ID, jun_poss in by_gene_ID.iteritems():
#                 jun_poss.sort(reverse=(strand == '-'))
#                 i_from_start = 1
#                 i_from_end = len(jun_poss)
#                 for jun_pos in jun_poss:
#                     assert tmpd.setdefault(jun_pos,
#                                            (i_from_start, i_from_end)) == (
#                            i_from_start, i_from_end)
#                     i_from_start += 1
#                     i_from_end -= 1
#             self.jun_pos2gene_order[(chrome, strand)] = tmpd
#
#         # order junctions records by increasing junction start position
#         self._bisect_juncs_sta = {}
#         for (chrome, strand), jun_recs in self.junctions.iteritems():
#             jun_recs.sort()
#             juns_start = [x[0] for x in jun_recs]
#             self._bisect_juncs_sta[(chrome, strand)] = juns_start
#
#         # check that junction section ends here
#         r = f.readline().rstrip('\n\r')
#         assert (r == '')
#
#         ## junction distribution
#         self.distrib = {}
#         poss = f.readline().rstrip('\n\r').split('\t')
#         vals = f.readline().rstrip('\n\r').split('\t')
#         assert (len(poss) == len(vals))
#         if poss == ['']:
#             return
#         for p, v in zip(poss, vals):
#             p = int(p)
#             if self.maxU <= p <= self.maxD:
#                 self.distrib[p] = int(v)
#
#     def set_current_chrome_strand(self, chrome, strand):
#         anti_strand = ['+', '-'][strand == '+']
#         self._cur_strand = strand
#         self._cur_strand_fact = [-1, 1][strand == '+']
#         self._cur_juncs_same = self.junctions.get((chrome, strand), [])
#         self._cur_juncs_anti = self.junctions.get((chrome, anti_strand), [])
#         self._cur_bisect_juncs_sta_same = self._bisect_juncs_sta.get(
#             (chrome, strand), [])
#         self._cur_bisect_juncs_sta_anti = self._bisect_juncs_sta.get(
#             (chrome, anti_strand), [])
#
#     def get_closest(self, position):
#         """" Return distance relative to strand on which hit resides (p_start)
#              return junction_id, orient, d"""
#         jun_rec_same = None
#         if self._cur_bisect_juncs_sta_same:
#             i = bisect.bisect(self._cur_bisect_juncs_sta_same, position) - 1
#             if i >= 0:
#                 hit_jun_rec = self._cur_juncs_same[i]
#                 hit_jun_start, hit_jun_stop, _, _, hit_jun_pos, _, _, _, _ = hit_jun_rec
#                 if position <= hit_jun_stop:
#                     tmp_d = (position - hit_jun_pos) * self._cur_strand_fact
#                     if self.maxU <= tmp_d <= self.maxD:
#                         jun_rec_same = hit_jun_rec
#                         d_same = tmp_d
#         if jun_rec_same is None:
#             d_same = None
#
#         jun_rec_anti = None
#         if self._cur_bisect_juncs_sta_anti:
#             i = bisect.bisect(self._cur_bisect_juncs_sta_anti, position) - 1
#             if i >= 0:
#                 hit_jun_rec = self._cur_juncs_anti[i]
#                 hit_jun_start, hit_jun_stop, _, _, hit_jun_pos, _, _, _, _ = hit_jun_rec
#                 if position <= hit_jun_stop:
#                     tmp_d = (hit_jun_pos - position) * self._cur_strand_fact
#                     if self.maxU <= tmp_d <= self.maxD:
#                         jun_rec_anti = hit_jun_rec
#                         d_anti = tmp_d
#         if jun_rec_anti is None:
#             d_anti = None
#
#         if d_same is not None and d_anti is not None:
#             if abs(d_same) <= abs(d_anti):
#                 return jun_rec_same, 'same', d_same
#             else:
#                 return jun_rec_anti, 'anti', d_anti
#         elif d_same is not None:
#             return jun_rec_same, 'same', d_same
#         elif d_anti is not None:
#             return jun_rec_anti, 'anti', d_anti
#         else:
#             return None, None, None
#
#
# def generate_maps(bed_d, mapped_to, annot_ver, maxU_within, maxD_within,
#                   maxU_flanking, maxD_flanking, selected_chromes,
#                   analysis_id=None):
#     junc_to_load = [_fname_RNAmap(mapped_to, annot_ver, mname) for
#                     mname, _, _, _, _, _, _, _ in maps_to_gen]
#     #    land_to_load = glob.glob(os.path.join(iCount.storage_root, "genomes/annotation/RNAmap_la
# ndmarks/%s_%s_RNAmap_*.tab" % (mapped_to, annot_ver)))
#     land_to_load = []
#     land_to_load.sort()
#     maps_to_load = junc_to_load + land_to_load
#     # RNAmap object,
#     # map_xlink_same, map_xlink_anti, map_value_same, map_value_anti, map_junctions_same, map_jun
# ctions_anti,
#     # junction_xlinks_same, junction_xlinks_anti, junction_value_same, junction_value_anti,
#     # list_of_xlinks_and_distances_to_junction
#     rna_maps = [
#         (RNAmap(mn, maxU_within, maxD_within, maxU_flanking, maxD_flanking),
#          {}, {}, {}, {}, {}, {},
#          {}, {}, {}, {},
#          []) for mn in maps_to_load]
#
#     # start counting
#     xlinks_total = 0
#     count_total = 0
#     xlinks_used = {}  # for each map
#     count_used = {}  # for each map
#     xlinks_used_in_any = 0  # used in any map
#     count_used_in_any = 0  # used in any map
#     tot = len(bed_d)
#     prev_prog = 5
#     prog_i = 0
#     for (chrome, strand), tmpd in bed_d.iteritems():
#         if selected_chromes <> 'whole genome' and chrome not in selected_chromes:
#             # skip chromosome not selected for analysis
#             continue
#         # store "pointers" to landmarks map for current chrome and strand so
#         # that we do not have to look for it at each position
#         [x[0].set_current_chrome_strand(chrome, strand) for x in rna_maps]
#         for (p_start, p_end), value in tmpd.iteritems():
#             xlinks_total += 1
#             count_total += value
#
#             used = False
#             for (rna_map, map_xlink_same, map_xlink_anti,
#                  map_value_same, map_value_anti, map_junctions_same,
#                  map_junctions_anti,
#                  junction_xlinks_same, junction_xlinks_anti,
#                  junction_value_same, junction_value_anti,
#                  used_junction_xlinks_list) in rna_maps:
#                 # find closest junction
#                 junction_id, orient, d = rna_map.get_closest(p_start)
#
#                 if junction_id is None:
#                     assert (orient is None)
#                     assert (d is None)
#                     continue
#
#                 # redefine junction_id to be used subsequently
#                 jun_start, jun_stop, jun_chrome, jun_strand, jun_pos, spanU, spanD, associated_ge
# ne_name, associated_gene_ID = junction_id
#                 junction_id = (jun_chrome, jun_strand, jun_pos, spanU, spanD,
#                                associated_gene_name, associated_gene_ID)
#
#                 used = True
#                 if orient == 'same':
#                     map_xlink_same[d] = map_xlink_same.get(d, 0) + 1
#                     map_value_same[d] = map_value_same.get(d, 0) + value
#                     map_junctions_same.setdefault(d, set()).add(junction_id)
#                     tmpd = junction_xlinks_same.setdefault(junction_id, {})
#                     tmpd[d] = tmpd.get(d, 0) + 1
#                     tmpd = junction_value_same.setdefault(junction_id, {})
#                     tmpd[d] = tmpd.get(d, 0) + value
#                 else:
#                     map_xlink_anti[d] = map_xlink_anti.get(d, 0) + 1
#                     map_value_anti[d] = map_value_anti.get(d, 0) + value
#                     map_junctions_anti.setdefault(d, set()).add(junction_id)
#                     tmpd = junction_xlinks_anti.setdefault(junction_id, {})
#                     tmpd[d] = tmpd.get(d, 0) + 1
#                     tmpd = junction_value_anti.setdefault(junction_id, {})
#                     tmpd[d] = tmpd.get(d, 0) + value
#
#                 used_junction_xlinks_list.append((chrome, p_start, p_end,
#                                                   strand, value, junction_id,
#                                                   orient, d))
#                 xlinks_used[rna_map.name] = xlinks_used.get(rna_map.name,
#                                                             0) + 1
#                 count_used[rna_map.name] = count_used.get(rna_map.name,
#                                                           0) + value
#
#             if used:
#                 xlinks_used_in_any += 1
#                 count_used_in_any += value
#
#         if analysis_id is not None:
#             prog_i += 1
#             prog = 5 + (100000 * prog_i / tot) / 1190
#             if prog > prev_prog:
#                 db.analysis.set_status(analysis_id, 'processing %s%%' % prog)
#                 prev_prog = prog
#
#     return rna_maps, xlinks_total, xlinks_used, xlinks_used_in_any, count_total, count_used, coun
# t_used_in_any
#
#
# def smooth(data, hws, maxU=None, maxD=None):
#     def w_val(d, v, f, variance):
#         return v * f * math.exp(-0.5 * d * d / variance)
#
#     s = hws / 2.0  # /4.0
#     variance = s * s
#     f = 1.0 / math.sqrt(2.0 * math.pi * variance)
#
#     all_pos = data.keys()
#     if maxU is not None:
#         stop_at_maxU = maxU + hws
#         all_pos = [p for p in all_pos if p >= stop_at_maxU]
#     if maxD is not None:
#         stop_at_maxD = maxD - hws
#         all_pos = [p for p in all_pos if p <= stop_at_maxD]
#
#     ndata = {}
#     for p in all_pos:
#         nv = [w_val(p - p2, data.get(p2, 0.0), f, variance) for p2 in
#               range(p - hws, p + hws + 1)]
#         ndata[p] = sum(nv) / float(len(nv))
#
#     # rescale so that area under curve (absolute counts) is same as in original data
#     sum_data = sum([data.get(p, 0.0) for p in all_pos])
#     sum_ndata = sum(ndata.values())
#
#     if sum_ndata <= 0.0:
#         return ndata
#
#     f = 1.0 * sum_data / sum_ndata
#     n_ndata = {}
#     for p, v in ndata.iteritems():
#         n_ndata[p] = f * v
#     return n_ndata
#
#
# def normalize(rna_map_distrib, data, additional_factor=1.0):
#     ret = {}
#     for p, v in data.iteritems():
#         nf = rna_map_distrib.get(p, 0)
#         if nf:
#             ret[p] = additional_factor * float(v) / nf
#         else:
#             ret[p] = 0
#     return ret
#
#
# def draw_map(plt, data_same, data_anti, title, maxU, maxD, color_same,
#              color_anti):
#     plt.title(title)
#
#     # same strand
#     vals = [(p, v) for p, v in data_same.iteritems()]
#     if vals:
#         vals.sort()
#         posx_same, xlinks_same = zip(*vals)
#         poss_same = [p - 0.4 for p in posx_same]
#         plt.bar(poss_same, xlinks_same, width=0.8, color=color_same,
#                 edgecolor=color_same)
#         t_y_max = max(xlinks_same)
#     else:
#         t_y_max = 0
#     # anti strand (drawn as negative values, and in red)
#     vals = [(p, v) for p, v in data_anti.iteritems()]
#     if vals:
#         vals.sort()
#         posx_anti, xlinks_anti = zip(*vals)
#         xlinks_anti = [-v for v in xlinks_anti]
#         poss_anti = [p - 0.35 for p in posx_anti]
#         plt.bar(poss_anti, xlinks_anti, width=0.7, color=color_anti,
#                 edgecolor=color_anti)
#         t_y_min = min(xlinks_anti)
#     else:
#         t_y_min = 0
#
#     plt.plot([maxU, maxD], [0, 0], 'k', lw=0.1)
#     plt.plot([0, 0], [t_y_min * 1.01, t_y_max * 1.01], 'k', lw=0.2)
#     plt.xlim((maxU - 1, maxD + 1))
#     plt.ylim((t_y_min * 1.01, t_y_max * 1.01))
#
#
# def render_maps_text_and_html(bed_res_type, rna_maps, mapped_to,
#                               smoothing_factor, xlinks_total, xlinks_used,
#                               xlinks_used_in_any, count_total, count_used,
#                               count_used_in_any, out_fname_pref,
#                               color_same="green", color_anti="red", heading=[],
#                               analysis_id=None):
#     ret_files = []
#
#     ### overview stats ###################
#     # create output folder if does not exist yet
#     iCount.files.osutils.make_dir(os.path.split(out_fname_pref)[0])
#
#     # _maps.html
#     out_fname_html = os.path.join(os.path.split(out_fname_pref)[0],
#                                   "0main_" + os.path.split(out_fname_pref)[
#                                       1] + "_maps.html")
#     ret_files.append(('html-maps', out_fname_html))
#     fout_html = iCount.files.nc_open(out_fname_html, "wt")
#     # _stats.tab.gz
#     out_fname_stats = out_fname_pref + "_stats.tab.gz"
#     ret_files.append(('tab-stats', out_fname_stats))
#     fout_stats = iCount.files.nc_open(out_fname_stats, "wt")
#     # _maps.tab.gz
#     out_fname_maps = out_fname_pref + "_maps.tab.gz"
#     ret_files.append(('tab-maps', out_fname_maps))
#     fout_maps = iCount.files.nc_open(out_fname_maps, "wt")
#
#     fout_html.write('<html><head>\n')
#     fout_html.write('<style type="text/css">\n')
#     fout_html.write('td {font-size: 10pt;}\n')
#     fout_html.write('th {font-size: 10pt;}\n')
#     fout_html.write('</style>\n')
#
#     fout_html.write('<b>Parameters</b><br/>\n')
#     link_name = os.path.basename(out_fname_pref)
#     for h_p, h_v in heading:
#         fout_html.write('%s: %s<br/>\n' % (h_p, h_v))
#         fout_stats.write('%s\t%s\n' % (h_p, h_v))
#     fout_html.write('<br/>\n')
#     fout_stats.write('\n')
#     # distribution of x-links and vals
#     fout_html.write(
#         '<b>Distribution of x-link sites and counts in RNAmaps.</b><br/>\n')
#     fout_html.write('[x-link sites] total: %s<br/>\n' % xlinks_total)
#     fout_stats.write('[x-link sites] total\t%s\n' % xlinks_total)
#     if xlinks_total == 0.0: xlinks_total = 1.0
#     xlinks_used_in_any_per = 100.0 * xlinks_used_in_any / xlinks_total
#     fout_html.write('[x-link sites] used in any map: %g (%0.2f%%)<br/>\n' % (
#     xlinks_used_in_any, xlinks_used_in_any_per))
#     fout_stats.write('[x-link sites] used in any map\t%g\t%0.2f%%\n' % (
#     xlinks_used_in_any, xlinks_used_in_any_per))
#     if xlinks_used_in_any == 0.0: xlinks_used_in_any = 1.0
#
#     fout_html.write('[%s] total: %s<br/>\n' % (bed_res_type, count_total))
#     fout_stats.write('[%s] total\t%s\n' % (bed_res_type, count_total))
#     if count_total == 0.0: count_total = 1.0
#     count_used_in_any_per = 100.0 * count_used_in_any / count_total
#     fout_html.write('[%s] used in any map: %g (%0.2f%%)<br/>\n' % (
#     bed_res_type, count_used_in_any, count_used_in_any_per))
#     fout_stats.write('[%s] used in any map\t%g\t%0.2f%%\n' % (
#     bed_res_type, count_used_in_any, count_used_in_any_per))
#     if count_used_in_any == 0.0: count_used_in_any = 1.0
#     fout_stats.write('\n')
#
#     # normalization by data complexity factors
#     norm_complexity_xlinks_total = 1000000000.0 / xlinks_total
#     norm_complexity_count_total = 1000000000.0 / count_total
#
#     fout_maps.write(
#         "[x-link sites] complexity normalization factor (10^9 / total [x-link sites]): %g\n" % (
#         norm_complexity_xlinks_total))
#     fout_maps.write(
#         "[%s] complexity normalization factor (= 10^9 / total [%s]): %g\n\n" % (
#         bed_res_type, bed_res_type, norm_complexity_count_total))
#
#     fout_html.write('<table border=1>')
#     fout_html.write('<tr><th>map</th><th>span<br/>(upstream; downstream)</th>')
#     fout_html.write(
#         '<th>[x-link sites]<br/>used in map</th><th>[x-link sites]<br/>% of total</th><th>[x-link
#  sites]<br/>% of used in any map</th>')
#     fout_html.write(
#         '<th>[%s]<br/>used in map</th><th>[%s]<br/>%% of total</th><th>[%s]<br/>%% of used in any
#  map</th></tr>\n' % (
#         bed_res_type, bed_res_type, bed_res_type))
#
#     fout_stats.write('map\tupstream\tdownstream\t')
#     fout_stats.write(
#         '[x-link sites] used in map\t[x-link sites] % of total\t[x-link sites] % of used in any m
# ap\t')
#     fout_stats.write(
#         '[%s] used in map\t[%s] %% of total\t[%s] %% of used in any map\n' % (
#         bed_res_type, bed_res_type, bed_res_type))
#
#     for x in rna_maps:
#         name = x[0].name
#         cn_xlink_used = xlinks_used.get(name, 0.0)
#         cn_xlink_per_total = 100.0 * cn_xlink_used / xlinks_total
#         cn_xlink_per_used_in_any = 100.0 * cn_xlink_used / xlinks_used_in_any
#         cn_count_used = count_used.get(name, 0.0)
#         cn_count_per_total = 100.0 * cn_count_used / count_total
#         cn_count_per_used_in_any = 100.0 * cn_count_used / count_used_in_any
#         fout_html.write(
#             '<tr><td>%s</td><td>%s; %s</td>' % (name, x[0].maxU, x[0].maxD))
#         fout_html.write('<td>%g</td><td>%0.2f%%</td><td>%0.2f%%</td>' % (
#         cn_xlink_used, cn_xlink_per_total, cn_xlink_per_used_in_any))
#         fout_html.write(
#             '<td>%g</td><td>%0.2f%%</td><td>%0.2f%%</td></tr>\n' % (
#             cn_count_used, cn_count_per_total, cn_count_per_used_in_any))
#
#         fout_stats.write('%s\t%s\t%s\t' % (name, x[0].maxU, x[0].maxD))
#         fout_stats.write('%g\t%0.2f%%\t%0.2f%%\t' % (
#         cn_xlink_used, cn_xlink_per_total, cn_xlink_per_used_in_any))
#         fout_stats.write('%g\t%0.2f%%\t%0.2f%%\n' % (
#         cn_count_used, cn_count_per_total, cn_count_per_used_in_any))
#     fout_html.write('</table>\n')
#     fout_html.write(
#         '<a href="%s">Download statistics in tab-delimited format.</a><br/>\n' % os.path.basename
# (
#             out_fname_stats))
#     fout_html.write(
#         '<a href="%s">Download graphs in tab-delimited format.</a><br/>\n' % os.path.basename(
#             out_fname_maps))
#     fout_html.write('<hl>\n')
#     fout_stats.write('\n')
#
#     num_columns = 3
#
#     ### individual maps ##################
#     if smoothing_factor > 0:
#         plt_num_rows = 4
#         plt.figure(1, figsize=(num_columns * 4, 9.72))
#         plt.subplots_adjust(left=0.05, right=0.99, bottom=0.035, top=0.964,
#                             wspace=0.22, hspace=0.22)
#     else:
#         plt_num_rows = 2
#         plt.figure(1, figsize=(num_columns * 4, 5))
#         plt.subplots_adjust(left=0.05, right=0.99, bottom=0.07, top=0.93,
#                             wspace=0.22, hspace=0.22)
#     tot = len(rna_maps)
#     prev_prog = 90
#     prog_i = 0
#
#     for (rna_map,
#          map_xlink_same, map_xlink_anti, map_value_same, map_value_anti,
#          map_junctions_same, map_junctions_anti,
#          junction_xlinks_same, junction_xlinks_anti, junction_value_same,
#          junction_value_anti,
#          used_junction_xlinks_list) in rna_maps:
#         #
#         plt.clf()
#         xlabel = "position relative to %s" % rna_map.kind
#
#         # x-links - raw
#         plt.subplot(plt_num_rows, num_columns,
#                     1 + num_columns * 0 + 0)  # (0,0)
#         draw_map(plt, map_xlink_same, map_xlink_anti,
#                  '[x-link sites] at position', rna_map.maxU, rna_map.maxD,
#                  color_same, color_anti)
#         # x-links normalized by junction distribution and data complexity
#         plt.subplot(plt_num_rows, num_columns,
#                     1 + num_columns * 0 + 2)  # 3) # (0,3)
#         complex_norm_map_xlink_same = normalize(rna_map.distrib,
#                                                 map_xlink_same,
#                                                 norm_complexity_xlinks_total)
#         complex_norm_map_xlink_anti = normalize(rna_map.distrib,
#                                                 map_xlink_anti,
#                                                 norm_complexity_xlinks_total)
#         draw_map(plt, complex_norm_map_xlink_same, complex_norm_map_xlink_anti,
#                  '10^9 * [x-link sites] / total [x-link sites]\nnormalized by all %ss spanning po
# sition' % (
#                  rna_map.kind), rna_map.maxU, rna_map.maxD, color_same,
#                  color_anti)
#
#         # values - raw
#         plt.subplot(plt_num_rows, num_columns,
#                     1 + num_columns * 1 + 0)  # (1,0)
#         draw_map(plt, map_value_same, map_value_anti,
#                  "[%s] at position" % (bed_res_type), rna_map.maxU,
#                  rna_map.maxD, color_same, color_anti)
#         if smoothing_factor == 0:
#             plt.xlabel(xlabel)
#         # values - normalized by junction distribution and data complexity
#         plt.subplot(plt_num_rows, num_columns,
#                     1 + num_columns * 1 + 2)  # 3) # (1,3)
#         complex_norm_map_value_same = normalize(rna_map.distrib,
#                                                 map_value_same,
#                                                 norm_complexity_count_total)
#         complex_norm_map_value_anti = normalize(rna_map.distrib,
#                                                 map_value_anti,
#                                                 norm_complexity_count_total)
#         draw_map(plt, complex_norm_map_value_same, complex_norm_map_value_anti,
#                  "10^9 * [%s] / total [%s]\nnormalized by all %ss spanning position" % (
#                  bed_res_type, bed_res_type, rna_map.kind), rna_map.maxU,
#                  rna_map.maxD, color_same, color_anti)
#         if smoothing_factor == 0:
#             plt.xlabel(xlabel)
#
#         # junctions (landmarks)
#         plt.subplot(plt_num_rows, num_columns,
#                     1 + num_columns * 0 + 1)  # (0, 1)
#         plt.title("%s %ss\ndistribution of all %ss spanning position" % (
#         rna_map.name, rna_map.kind, rna_map.kind))
#         vals = [(p, v) for p, v in rna_map.distrib.iteritems() if
#                 rna_map.maxU <= p <= rna_map.maxD]
#         if vals:
#             vals.sort()
#             posv, junctions_same = zip(*vals)
#             l_same = plt.bar(posv, junctions_same, width=0.8, color='gray',
#                              edgecolor='gray')
#             t_y_max = max(junctions_same)
#         else:
#             t_y_max = 0
#
#         t_y_min = 0
#         plt.plot([0, 0], [t_y_min * 1.01, t_y_max * 1.01], 'k', lw=0.2)
#
#         leg_lines = [
#             plt.plot([rna_map.maxU - 100], [0], color_same),
#             plt.plot([rna_map.maxU - 100], [0], color_anti),
#         ]
#         plt.ylim((t_y_min * 1.01, t_y_max * 1.01))
#         plt.xlim((rna_map.maxU - 1, rna_map.maxD + 1))
#         plt.figlegend([x[0] for x in leg_lines], ('same sense', 'anti sense'),
#                       loc='lower center')
#
#         # junctions are same for both measures
#         jun_all_same = set()
#         for v in map_junctions_same.values():
#             jun_all_same.update(v)
#         jun_all_anti = set()
#         for v in map_junctions_anti.values():
#             jun_all_anti.update(v)
#         jun_all_both = jun_all_same | jun_all_anti
#
#         cn_jun_same = float(len(jun_all_same))
#         cn_jun_anti = float(len(jun_all_anti))
#         cn_jun_both = float(len(jun_all_both))
#         cn_jun_both_ref = max(1.0, cn_jun_both)
#
#         # junctions (landmarks) of only junctions that were used
#         used_junctions_distrib = {}
#         for jun_chrome, jun_strand, jun_pos, jun_up, jun_down, jun_associated_gene_name, jun_asso
# ciated_gene_ID in jun_all_both:
#             for p in range(jun_up, jun_down + 1):
#                 used_junctions_distrib[p] = used_junctions_distrib.get(p,
#                                                                        0) + 1
#
#         plt.subplot(plt_num_rows, num_columns,
#                     1 + num_columns * 1 + 1)  # (1, 1)
#         plt.title(
#             "distribution of used %ss spanning position" % (rna_map.kind))
#         plt.xlabel("position relative to %s" % rna_map.kind)
#         vals = [(p, v) for p, v in used_junctions_distrib.iteritems() if
#                 rna_map.maxU <= p <= rna_map.maxD]
#         if vals:
#             vals.sort()
#             posv, junctions_same = zip(*vals)
#             l_same = plt.bar(posv, junctions_same, width=0.8, color='gray',
#                              edgecolor='gray')
#             t_y_max = max(junctions_same)
#         else:
#             t_y_max = 0
#
#         t_y_min = 0
#         plt.plot([0, 0], [t_y_min * 1.01, t_y_max * 1.01], 'k', lw=0.2)
#         plt.ylim((t_y_min * 1.01, t_y_max * 1.01))
#         plt.xlim((rna_map.maxU - 1, rna_map.maxD + 1))
#
#         # smooth maps
#         if smoothing_factor > 0:
#             smooth_map_xlink_same = smooth(map_xlink_same, smoothing_factor,
#                                            rna_map.maxU, rna_map.maxD)
#             smooth_map_xlink_anti = smooth(map_xlink_anti, smoothing_factor,
#                                            rna_map.maxU, rna_map.maxD)
#             smooth_map_value_same = smooth(map_value_same, smoothing_factor,
#                                            rna_map.maxU, rna_map.maxD)
#             smooth_map_value_anti = smooth(map_value_anti, smoothing_factor,
#                                            rna_map.maxU, rna_map.maxD)
#
#             # smooth function
#             plt.subplot(plt_num_rows, num_columns,
#                         1 + num_columns * 2 + 1)  # (2,1)
#             #            sm_fun = smooth(dict([(p, 1.0) for p in range(-smoothing_factor*3, smoot
# hing_factor*3+1)]), smoothing_factor)
#             sm_fun = {}
#             for p in range(-smoothing_factor * 3, smoothing_factor * 3 + 1):
#                 sm_fun[p] = 0.0
#             sm_fun[0] = 1.0
#             sm_fun = smooth(sm_fun, smoothing_factor)
#             draw_map(plt, sm_fun, {},
#                      'Smoothing function (half-width = %s)' % (
#                      smoothing_factor), -smoothing_factor - 1,
#                      smoothing_factor + 1, 'blue', 'white')
#             plt.xlabel('position relative to x-link')
#             plt.yticks([])
#             plt.ylabel('smoothing factor')
#
#             # x-links - raw
#             plt.subplot(plt_num_rows, num_columns,
#                         1 + num_columns * 2 + 0)  # (2,0)
#             draw_map(plt, smooth_map_xlink_same, smooth_map_xlink_anti,
#                      'smoothed [x-link sites] at position', rna_map.maxU,
#                      rna_map.maxD, color_same, color_anti)
#             # x-links normalized by junction distribution and data complexity
#             plt.subplot(plt_num_rows, num_columns,
#                         1 + num_columns * 2 + 2)  # 3) # (2,3)
#             complex_norm_smooth_map_xlink_same = normalize(rna_map.distrib,
#                                                            smooth_map_xlink_same,
#                                                            norm_complexity_xlinks_total)
#             complex_norm_smooth_map_xlink_anti = normalize(rna_map.distrib,
#                                                            smooth_map_xlink_anti,
#                                                            norm_complexity_xlinks_total)
#             draw_map(plt, complex_norm_smooth_map_xlink_same,
#                      complex_norm_smooth_map_xlink_anti,
#                      'smoothed 10^9 * [x-link sites] / total [x-link sites]\nnormalized by all %s
# s spanning position' % (
#                      rna_map.kind), rna_map.maxU, rna_map.maxD, color_same,
#                      color_anti)
#
#             # values - raw
#             plt.subplot(plt_num_rows, num_columns,
#                         1 + num_columns * 3 + 0)  # (3,0)
#             draw_map(plt, smooth_map_value_same, smooth_map_value_anti,
#                      "smoothed [%s] at position" % (bed_res_type),
#                      rna_map.maxU, rna_map.maxD, color_same, color_anti)
#             if smoothing_factor > 0:
#                 plt.xlabel(xlabel)
#             # values - normalized by junction distribution and data complexity
#             plt.subplot(plt_num_rows, num_columns,
#                         1 + num_columns * 3 + 2)  # 3) # (3,3)
#             complex_norm_smooth_map_value_same = normalize(rna_map.distrib,
#                                                            smooth_map_value_same,
#                                                            norm_complexity_count_total)
#             complex_norm_smooth_map_value_anti = normalize(rna_map.distrib,
#                                                            smooth_map_value_anti,
#                                                            norm_complexity_count_total)
#             draw_map(plt, complex_norm_smooth_map_value_same,
#                      complex_norm_smooth_map_value_anti,
#                      "smoothed 10^9 * [%s] / total [%s]\nnormalized by all %ss spanning position"
#  % (
#                      bed_res_type, bed_res_type, rna_map.kind), rna_map.maxU,
#                      rna_map.maxD, color_same, color_anti)
#             if smoothing_factor > 0:
#                 plt.xlabel(xlabel)
#
#         ## save maps to file in tabular format
#         fout_maps.write('RNAmap\t%s\n' % rna_map.name)
#         all_poss = set(rna_map.distrib.keys()) | set(
#             map_xlink_same.keys()) | set(map_xlink_anti.keys()) | set(
#             map_value_same.keys()) | set(map_value_anti.keys())
#
#         all_poss = [p for p in sorted(all_poss) if
#                     rna_map.maxU <= p <= rna_map.maxD]
#         fout_maps.write(
#             "position\t%s\n" % "\t".join([str(x) for x in all_poss]))
#
#         vals = ["%g" % rna_map.distrib.get(p, 0.0) for p in all_poss]
#         fout_maps.write("junction distribution\t%s\n" % "\t".join(vals))
#
#         vals = ["%g" % used_junctions_distrib.get(p, 0.0) for p in all_poss]
#         fout_maps.write("used junction distribution\t%s\n" % "\t".join(vals))
#
#         vals = ["%g" % (norm_complexity_xlinks_total / rna_map.distrib.get(p,
#                                                                            0.0) if p in rna_map.d
# istrib else 0.0)
#                 for p in all_poss]
#         fout_maps.write(
#             "[x-link sites] complexity normalization factor\t%s\n" % "\t".join(
#                 vals))
#
#         vals = ["%g" % (norm_complexity_count_total / rna_map.distrib.get(p,
#                                                                           0.0) if p in rna_map.di
# strib else 0.0)
#                 for p in all_poss]
#         fout_maps.write("[%s] complexity normalization factor\t%s\n" % (
#         bed_res_type, "\t".join(vals)))
#
#         vals = ["%g" % map_xlink_same.get(p, 0.0) for p in all_poss]
#         fout_maps.write("[x-link sites] (same)\t%s\n" % "\t".join(vals))
#         vals = ["%g" % map_xlink_anti.get(p, 0.0) for p in all_poss]
#         fout_maps.write("[x-link sites] (anti)\t%s\n" % "\t".join(vals))
#
#         vals = ["%g" % map_value_same.get(p, 0.0) for p in all_poss]
#         fout_maps.write("[%s] (same)\t%s\n" % (bed_res_type, "\t".join(vals)))
#         vals = ["%g" % map_value_anti.get(p, 0.0) for p in all_poss]
#         fout_maps.write("[%s] (anti)\t%s\n" % (bed_res_type, "\t".join(vals)))
#
#         vals = ["%g" % complex_norm_map_xlink_same.get(p, 0.0) for p in
#                 all_poss]
#         fout_maps.write(
#             "complexity normalized [x-link sites] (same)\t%s\n" % "\t".join(
#                 vals))
#         vals = ["%g" % complex_norm_map_xlink_anti.get(p, 0.0) for p in
#                 all_poss]
#         fout_maps.write(
#             "complexity normalized [x-link sites] (anti)\t%s\n" % "\t".join(
#                 vals))
#
#         vals = ["%g" % complex_norm_map_value_same.get(p, 0.0) for p in
#                 all_poss]
#         fout_maps.write("complexity normalized [%s] (same)\t%s\n" % (
#         bed_res_type, "\t".join(vals)))
#         vals = ["%g" % complex_norm_map_value_anti.get(p, 0.0) for p in
#                 all_poss]
#         fout_maps.write("complexity normalized [%s] (anti)\t%s\n" % (
#         bed_res_type, "\t".join(vals)))
#
#         if smoothing_factor > 0:
#             vals = ["%g" % smooth_map_xlink_same.get(p, 0.0) for p in all_poss]
#             fout_maps.write(
#                 "smoothed [x-link sites] (same)\t%s\n" % "\t".join(vals))
#             vals = ["%g" % smooth_map_xlink_anti.get(p, 0.0) for p in all_poss]
#             fout_maps.write(
#                 "smoothed [x-link sites] (anti)\t%s\n" % "\t".join(vals))
#
#             vals = ["%g" % smooth_map_value_same.get(p, 0.0) for p in all_poss]
#             fout_maps.write(
#                 "smoothed [%s] (same)\t%s\n" % (bed_res_type, "\t".join(vals)))
#             vals = ["%g" % smooth_map_value_anti.get(p, 0.0) for p in all_poss]
#             fout_maps.write(
#                 "smoothed [%s] (anti)\t%s\n" % (bed_res_type, "\t".join(vals)))
#
#             vals = ["%g" % complex_norm_smooth_map_xlink_same.get(p, 0.0) for p
#                     in all_poss]
#             fout_maps.write(
#                 "complexity normalized smoothed [x-link sites] (same)\t%s\n" % "\t".join(
#                     vals))
#             vals = ["%g" % complex_norm_smooth_map_xlink_anti.get(p, 0.0) for p
#                     in all_poss]
#             fout_maps.write(
#                 "complexity normalized smoothed [x-link sites] (anti)\t%s\n" % "\t".join(
#                     vals))
#
#             vals = ["%g" % complex_norm_smooth_map_value_same.get(p, 0.0) for p
#                     in all_poss]
#             fout_maps.write(
#                 "complexity normalized smoothed [%s] (same)\t%s\n" % (
#                 bed_res_type, "\t".join(vals)))
#             vals = ["%g" % complex_norm_smooth_map_value_anti.get(p, 0.0) for p
#                     in all_poss]
#             fout_maps.write(
#                 "complexity normalized smoothed [%s] (anti)\t%s\n" % (
#                 bed_res_type, "\t".join(vals)))
#         fout_maps.write("\n")
#
#         # save
#         plt.savefig("%s_%s.png" % (out_fname_pref, rna_map.name))
#         os.chmod("%s_%s.png" % (out_fname_pref, rna_map.name),
#                  stat.S_IRGRP | stat.S_IRUSR | stat.S_IROTH | stat.S_IWUSR)
#         plt.savefig("%s_%s.eps" % (out_fname_pref, rna_map.name))
#         os.chmod("%s_%s.eps" % (out_fname_pref, rna_map.name),
#                  stat.S_IRGRP | stat.S_IRUSR | stat.S_IROTH | stat.S_IWUSR)
#
#         link_name = os.path.basename("%s_%s" % (out_fname_pref, rna_map.name))
#         fout_html.write('<h2>RNAmap %s</h2>\n' % (
#         '"%s"' % rna_map.name if " " in rna_map.name else rna_map.name))
#         fout_html.write(
#             "<b>%s</b> is at position zero.<br/>\n" % (rna_map.kind))
#         fout_html.write("upstream is <b>%s</b><br/>\n" % (rna_map.descU))
#         fout_html.write("downstream is <b>%s</b><br/>\n" % (rna_map.descD))
#         fout_html.write(
#             '<a href="%s.eps"><img border="1" src="%s.png"></a>\n' % (
#             link_name, link_name))
#         fout_html.write('<a href="%s.eps">[eps]</a><br/><br/>\n' % (link_name))
#
#         fout_html.write(
#             '<small><b>x-links and %ss distributions</b></small><br/>\n' % (
#             rna_map.kind))
#
#         fout_html.write('<table border=1>\n')
#
#         # x-links and counts (output formatting)
#         frmts = ["%s", "%g", "%0.2f%%", "%0.2f%%", "%0.2f%%", "%g", "%g"]
#         frmt_html = "</td><td>".join(frmts)
#         frmt_txt = "\t".join(frmts)
#
#         # x-links
#         fout_html.write('<tr><th>strand</th>\n')
#         fout_html.write(
#             '<th>[x-link sites]<br/>used in map</th><th>[x-link sites]<br/>% of total</th><th>[x-
# link sites]<br/>% of used in any map</th><th>[x-link sites]<br/>% of used in map</th>')
#         fout_html.write(
#             '<th>unique [%ss]<br/>in map</th><th>[x-link sites]<br/>per unique [%s]</th></tr>\n'
# % (
#             rna_map.kind, rna_map.kind))
#         fout_stats.write('\nRNAmap\t%s\n' % rna_map.name)
#         fout_stats.write('strand\t')
#         fout_stats.write(
#             '[x-link sites] used in map\t[x-link sites] % of total\t[x-link sites] % of used in a
# ny map\t[x-link sites] % of used in map\t')
#         fout_stats.write(
#             'unique [%ss] in map\t[x-link sites] per unique [%s]\n' % (
#             rna_map.kind, rna_map.kind))
#
#         cn_used_xlink = xlinks_used.get(rna_map.name, 0)
#         if cn_used_xlink == 0.0: cn_used_xlink = 1.0
#         cn_used_in_any_xlink = xlinks_used_in_any
#         if cn_used_in_any_xlink == 0.0: cn_used_in_any_xlink = 1.0
#         cn_same_xlink = sum([0.0] + map_xlink_same.values())
#         cn_anti_xlink = sum(map_xlink_anti.values())
#         cn_both_xlink = cn_same_xlink + cn_anti_xlink
#
#         cn_xlink_per_jun_same = cn_same_xlink / cn_jun_both_ref
#         cn_xlink_per_jun_anti = cn_anti_xlink / cn_jun_both_ref
#         cn_xlink_per_jun_both = cn_both_xlink / cn_jun_both_ref
#
#         same_vals = (
#         'same', cn_same_xlink, 100.0 * cn_same_xlink / xlinks_total,
#         100.0 * cn_same_xlink / cn_used_in_any_xlink,
#         100.0 * cn_same_xlink / cn_used_xlink, cn_jun_same,
#         cn_xlink_per_jun_same)
#         anti_vals = (
#         'anti', cn_anti_xlink, 100.0 * cn_anti_xlink / xlinks_total,
#         100.0 * cn_anti_xlink / cn_used_in_any_xlink,
#         100.0 * cn_anti_xlink / cn_used_xlink, cn_jun_anti,
#         cn_xlink_per_jun_anti)
#         both_vals = (
#         'both', cn_both_xlink, 100.0 * cn_both_xlink / xlinks_total,
#         100.0 * cn_both_xlink / cn_used_in_any_xlink,
#         100.0 * cn_both_xlink / cn_used_xlink, cn_jun_both,
#         cn_xlink_per_jun_both)
#
#         fout_html.write("<tr><td>%s</td></tr>\n" % (frmt_html % same_vals))
#         fout_html.write("<tr><td>%s</td></tr>\n" % (frmt_html % anti_vals))
#         fout_html.write("<tr><td>%s</td></tr>\n" % (frmt_html % both_vals))
#         fout_stats.write("%s\n" % (frmt_txt % same_vals))
#         fout_stats.write("%s\n" % (frmt_txt % anti_vals))
#         fout_stats.write("%s\n" % (frmt_txt % both_vals))
#
#         # counts
#         fout_html.write('<tr><th>strand</th>\n')
#         fout_html.write(
#             '<th>[%s]<br/>used in map</th><th>[%s]<br/>%% of total</th><th>[%s]<br/>%% of used in
#  any map</th><th>[%s]<br/>%% of used in map</th>' % (
#             bed_res_type, bed_res_type, bed_res_type, bed_res_type))
#         fout_html.write(
#             '<th>unique [%ss]<br/>in map</th><th>[%s]<br/>per unique [%s]</th></tr>\n' % (
#             rna_map.kind, bed_res_type, rna_map.kind))
#
#         fout_stats.write('\nstrand\t')
#         fout_stats.write(
#             '[%s] used in map\t[%s] %% of total\t[%s] %% of used in any map\t[%s] %% of used in m
# ap\t' % (
#             bed_res_type, bed_res_type, bed_res_type, bed_res_type))
#         fout_stats.write('unique [%ss] in map\t[%s] per unique [%s]\n' % (
#         rna_map.kind, bed_res_type, rna_map.kind))
#
#         cn_used_value = count_used.get(rna_map.name, 0)
#         if cn_used_value == 0.0: cn_used_value = 1.0
#         cn_used_in_any_value = count_used_in_any
#         if cn_used_in_any_value == 0.0: cn_used_in_any_value = 1.0
#         cn_same_value = sum([0.0] + map_value_same.values())
#         cn_anti_value = sum([0.0] + map_value_anti.values())
#         cn_both_value = cn_same_value + cn_anti_value
#
#         cn_value_per_jun_same = cn_same_value / cn_jun_both_ref
#         cn_value_per_jun_anti = cn_anti_value / cn_jun_both_ref
#         cn_value_per_jun_both = cn_both_value / cn_jun_both_ref
#
#         same_vals = (
#         'same', cn_same_value, 100.0 * cn_same_value / count_total,
#         100.0 * cn_same_value / cn_used_in_any_value,
#         100.0 * cn_same_value / cn_used_value, cn_jun_same,
#         cn_value_per_jun_same)
#
#         anti_vals = (
#         'anti', cn_anti_value, 100.0 * cn_anti_value / count_total,
#         100.0 * cn_anti_value / cn_used_in_any_value,
#         100.0 * cn_anti_value / cn_used_value, cn_jun_anti,
#         cn_value_per_jun_anti)
#
#         both_vals = (
#         'both', cn_both_value, 100.0 * cn_both_value / count_total,
#         100.0 * cn_both_value / cn_used_in_any_value,
#         100.0 * cn_both_value / cn_used_value, cn_jun_both,
#         cn_value_per_jun_both)
#
#         fout_html.write("<tr><td>%s</td></tr>\n" % (frmt_html % same_vals))
#         fout_html.write("<tr><td>%s</td></tr>\n" % (frmt_html % anti_vals))
#         fout_html.write("<tr><td>%s</td></tr>\n" % (frmt_html % both_vals))
#         fout_stats.write("%s\n" % (frmt_txt % same_vals))
#         fout_stats.write("%s\n" % (frmt_txt % anti_vals))
#         fout_stats.write("%s\n" % (frmt_txt % both_vals))
#
#         fout_html.write('</table>\n')
#         fout_html.write('<br/>\n')
#
#         # render lists of top genes
#         def render_junction_sum_table(fout, jun_list, html=True):
#             if html:
#                 r = ["position of<br/>[%s]" % rna_map.kind,
#                      "span<br/>(upstream; downstream)",
#                      "associated gene_name",
#                      "associated gene_ID",
#                      "[x-link sites]<br/>(same; anti)",
#                      "[%s]<br/>(same; anti)" % bed_res_type,
#                      ]
#                 fout.write('<table border=1>')
#                 ostr = "<tr><th>" + "</th><th>".join(r) + "</th></tr>"
#             else:
#                 r = ["position of [%s]" % rna_map.kind,
#                      "span (upstream)",
#                      "span (downstream)",
#                      "associated gene_name",
#                      "associated gene_ID",
#                      "sum [x-link sites] (same)",
#                      "sum [x-link sites] (anti)",
#                      "sum [%s] (same)" % bed_res_type,
#                      "sum [%s] (anti)" % bed_res_type,
#                      ]
#                 ostr = "\t".join([x.replace('<br/>', '') for x in r])
#             fout.write("%s\n" % ostr)
#
#             for js, jun_id in jun_list:
#                 jun_chrome, jun_strand, jun_pos, jun_up, jun_down, jun_associated_gene_name, jun_
# associated_gene_ID = jun_id
#
#                 xl_same = sum(
#                     [0.0] + junction_xlinks_same.get(jun_id, {}).values())
#                 xl_anti = sum(
#                     [0.0] + junction_xlinks_anti.get(jun_id, {}).values())
#                 v_same = sum(
#                     [0.0] + junction_value_same.get(jun_id, {}).values())
#                 v_anti = sum(
#                     [0.0] + junction_value_anti.get(jun_id, {}).values())
#
#                 if html:
#                     r = ["%s:%s@%s" % (jun_chrome, jun_pos, jun_strand),
#                          "%s; %s" % (jun_up, jun_down),
#                          jun_associated_gene_name, jun_associated_gene_ID,
#                          "%g; %g" % (xl_same, xl_anti),
#                          "%g; %g" % (v_same, v_anti)]
#                     ostr = "<tr><td>" + "</td><td>".join(r) + "</td></tr>"
#                 else:
#                     r = ["%s:%s@%s" % (jun_chrome, jun_pos, jun_strand),
#                          "%s" % jun_up, "%s" % jun_down,
#                          jun_associated_gene_name, jun_associated_gene_ID,
#                          "%g" % xl_same, "%g" % xl_anti, "%g" % v_same,
#                          "%g" % v_anti]
#                     ostr = "\t".join(r)
#                 fout.write("%s\n" % ostr)
#             if html:
#                 fout.write("</table>\n")
#
#         take_top = 10
#         fout_html.write('<table>\n')
#         fout_html.write('<tr><td valign="top">\n')
#         fout_html.write(
#             "<b>Top %s %ss with largest x-link sites sum (both strands)</b><br/>\n" % (
#             take_top, rna_map.kind))
#         ks = set(junction_xlinks_same.keys()) | set(
#             junction_xlinks_anti.keys())
#         top_list = [(sum([0.0] + junction_xlinks_same.get(jun_id,
#                                                           {}).values() + junction_xlinks_anti.get
# (
#             jun_id, {}).values()), jun_id) for jun_id in ks]
#         top_list.sort(reverse=True)
#         top_list = top_list[:take_top]
#         top_list = top_list[:take_top]
#         render_junction_sum_table(fout_html, top_list)
#         fout_html.write('</td><td valign="top">\n')
#
#         fout_html.write(
#             "<b>Top %s %ss with largest [%s] sum (both strands)</b><br/>\n" % (
#             take_top, rna_map.kind, bed_res_type))
#         ks = set(junction_value_same.keys()) | set(junction_value_anti.keys())
#         top_list = [(sum([0.0] + junction_value_same.get(jun_id,
#                                                          {}).values() + junction_value_anti.get(
#             jun_id, {}).values()), jun_id) for jun_id in ks]
#         top_list.sort(reverse=True)
#         render_junction_sum_table(fout_html, top_list[
#                                              :take_top])  # need complete list to store it to fil
# e
#         fout_html.write("</td></tr>\n")
#         fout_html.write("</table>\n")
#
#         ### list of junctions in map ##########
#         out_fname_junct_list = out_fname_pref + "_sumjunct_%s.tab.gz" % (
#         rna_map.name)
#         fout_junct_list = iCount.files.nc_open(out_fname_junct_list, "wt")
#         ret_files.append(
#             ('tab-junctions_in_%s_map' % rna_map.name, out_fname_junct_list))
#         render_junction_sum_table(fout_junct_list, top_list, html=False)
#         fout_junct_list.close()
#
#         fout_html.write(
#             '<a href="%s">Download full junction summary list in tab-delimited format.</a><br/>\n
# ' % os.path.basename(
#                 out_fname_junct_list))
#
#         ### subset of x-links used for RNA-map, in two formats:
#         ### tab - with additional information on distance from junction, junction ID
#         ### bed - in bedGraph format, for viewing in UCSC, k-mer analysis, etc...
#         out_fname_junct_xlinks_list_tab = out_fname_pref + "_xlinks_annotated_%s.tab.gz" % (
#         rna_map.name)
#         out_fname_junct_xlinks_list_bed = out_fname_pref + "_xlinks_%s.bed.gz" % (
#         rna_map.name)
#         fout_junct_xlinks_list_tab = iCount.files.nc_open(
#             out_fname_junct_xlinks_list_tab, "wt")
#         fout_junct_xlinks_list_bed = iCount.files.nc_open(
#             out_fname_junct_xlinks_list_bed, "wt")
#         ret_files.append(('tab-junction_xlinks_in_%s_map' % rna_map.name,
#                           out_fname_junct_xlinks_list_tab))
#         ret_files.append(('bedGraph-junction_xlinks_in_%s_map' % rna_map.name,
#                           out_fname_junct_xlinks_list_bed))
#
#         used_junction_xlinks_list.sort()
#         tmp_bedGraph_data = []
#         fout_junct_xlinks_list_tab.write("%s\n" % "\t".join(
#             ["chromosome", "xlink_hit_position", "strand", "value",
#              "junction_ID", "associated gene_name", "associated gene_ID",
#              "xlink_orientation_relative_to_junction",
#              "xlink_distance_from_junction"]))
#         for (chrome, p_start, p_end, strand, value, junction_id, orient,
#              d) in used_junction_xlinks_list:
#             strand_value = "%s%s" % (strand, value)
#             value = "%g" % (value)
#             jun_chrome, jun_strand, jun_pos, jun_up, jun_down, jun_associated_gene_name, jun_asso
# ciated_gene_ID = junction_id
#             junction_id = "%s:%s@%s" % (jun_chrome, jun_pos, jun_strand)
#             olst = [chrome, str(p_start), strand, value, junction_id,
#                     jun_associated_gene_name, jun_associated_gene_ID, orient,
#                     str(d)]
#             ostr = "\t".join(olst)
#             fout_junct_xlinks_list_tab.write("%s\n" % ostr)
#             tmp_bedGraph_data.append((chrome, p_start, p_end, strand_value))
#         bg_description = "%s RNAmap (ID %s) x-links" % (
#         rna_map.name, analysis_id)
#         bg_name = os.path.basename(out_fname_junct_xlinks_list_bed)
#         BED.save_bedGraph_header(fout_junct_xlinks_list_bed, bg_name,
#                                  bg_description, mapped_to)
#         BED.save_bedGraph_records(fout_junct_xlinks_list_bed,
#                                   tmp_bedGraph_data)
#         fout_junct_xlinks_list_tab.close()
#         fout_junct_xlinks_list_bed.close()
#
#         fout_html.write(
#             '<a href="%s">Download full list of x-links used for RNAmap in tab-delimited format.<
# /a><br/>\n' % os.path.basename(
#                 out_fname_junct_xlinks_list_tab))
#
#         if analysis_id is not None:
#             prog_i += 1
#             prog = min(99, 90 + (100000 * prog_i / tot) / 10000)
#             if prog > prev_prog:
#                 db.analysis.set_status(analysis_id, 'processing %s%%' % prog)
#                 prev_prog = prog
#
#     plt.close(1)
#
#     fout_maps.close()
#     print
#     "data on individual RNAmaps stored to:", out_fname_maps
#     fout_stats.close()
#     print
#     "overview stats stored to:", out_fname_stats
#     fout_html.close()
#     print
#     "html report stored to:", out_fname_html
#     sys.stdout.flush()
#
#     return ret_files
#
#
# def read_results(fn):
#     res = {}
#     f = iCount.files.nc_open(fn, "rt")
#     # skip normalization information
#     r = f.readline()
#     r = f.readline()
#     r = f.readline()
#
#     r = f.readline()
#     while r:
#         r = r.rstrip('\n\r').split('\t')
#         if len(r) <> 2:
#             print >> sys.stderr, "Missing RNAmap name in file:", fn
#             print >> sys.stderr, "line:", r
#             return {}
#         if r[0] <> 'RNAmap':
#             print
#             sys.stderr, "Wrong RNAmap header format in file:", fn
#             print >> sys.stderr, "line:", r
#             return {}
#         map_name = r[1]
#         if map_name in res:
#             print
#             sys.stderr, "Duplicated RNAmap names in file:", fn
#             print >> sys.stderr, "line:", r
#             return {}
#         r = f.readline()
#         tmpd = {}
#         while r and r.rstrip('\n\r') <> '':
#             r = r.rstrip('\n\r').split('\t')
#             var_name = r[0]
#             vals = [float(x) for x in r[1:]]
#             tmpd[var_name] = vals
#             r = f.readline()
#         res[map_name] = tmpd
#         r = f.readline()
#     f.close()
#     return res
#
#
# def process_file(in_fname, mapped_to, annot_ver, out_fname_pref, maxU_within,
#                  maxD_within, maxU_flanking, maxD_flanking, smoothing_factor,
#                  selected_chromes, bed_res_type, analysis_id=None):
#     print
#     "*********************************"
#     print
#     time.ctime()
#     print
#     "generating RNAmaps for file:", in_fname
#     print
#     "saving to:", out_fname_pref
#     print
#     "parameters:"
#     print
#     "\tmax upstream within gene:", maxU_within
#     print
#     "\tmax downstream within gene:", maxD_within
#     print
#     "\tmax upstream gene flanking:", maxU_flanking
#     print
#     "\tmax downstream gene flanking:", maxD_flanking
#     print
#     "\thalf-width of Gaussian smoothing window:", smoothing_factor
#     print
#     "\tselected chromes:", selected_chromes
#     sys.stdout.flush()
#     if analysis_id is not None:
#         db.analysis.set_status(analysis_id, 'processing started')
#         prev_prog = 0
#
#     bed_d = BED.load_bedGraph(in_fname)
#     bed_res_type_desc = beds.desc_by_result_type.get(bed_res_type, None)
#     if bed_res_type_desc is not None:
#         bed_res_type = bed_res_type_desc[-1]
#     print
#     "input res_type: %s" % bed_res_type
#     sys.stdout.flush()
#
#     if analysis_id is not None:
#         db.analysis.set_status(analysis_id, 'processing 5%')
#         prev_prog = 5
#
#     rna_maps, xlinks_total, xlinks_used, xlinks_used_in_any, count_total, count_used, count_used_
# in_any = generate_maps(
#         bed_d, mapped_to, annot_ver, maxU_within, maxD_within, maxU_flanking,
#         maxD_flanking, selected_chromes, analysis_id=analysis_id)
#
#     print
#     "RNAmaps calculated"
#     sys.stdout.flush()
#
#     heading_params = [
#         ("input_file", os.path.basename(in_fname)),
#         ("input_res_type", bed_res_type),
#         ("mapped_to", mapped_to),
#         ("annot_ver", annot_ver),
#         ("maxU_within", maxU_within),
#         ("maxD_within", maxD_within),
#         ("maxU_flanking", maxU_flanking),
#         ("maxD_flanking", maxD_flanking),
#         ("smoothing", smoothing_factor),
#         ("selected_chromes", selected_chromes),
#     ]
#
#     files_generated = render_maps_text_and_html(bed_res_type, rna_maps,
#                                                 mapped_to, smoothing_factor,
#                                                 xlinks_total, xlinks_used,
#                                                 xlinks_used_in_any,
#                                                 count_total, count_used,
#                                                 count_used_in_any,
#                                                 out_fname_pref,
#                                                 heading=heading_params,
#                                                 analysis_id=analysis_id)
#
#     if analysis_id is not None:
#         files = [(ft, iCount.strip_storage_prefix(fn), db.filesize_nicefmt(fn))
#                  for ft, fn in files_generated]
#         db.analysis.set_output(analysis_id, files)
#         db.analysis.set_status(analysis_id, 'done')
#         db.analysis.set_last_update(analysis_id)
#     print
#     "DONE"
#     print
#     time.ctime()
#     print
#     "*********************************"
#     print
#     sys.stdout.flush()
#
#
# ### GENERATE JUNCTION RNAmaps from annotation
# def generate_RNAmap_defition_from_annotation():
#     for aname, annot in iCount.genomes.assembly.iteritems():
#         if aname <> 'hg19':
#             continue
#         for annot_ver, annot in iCount.genomes.assembly[
#             aname].annotation.iteritems():
#             if annot_ver == 'latest':
#                 continue
#
#             print
#             aname, annot_ver
#             sys.stdout.flush()
#             segmentation = annot.get_segmentation()
#             print
#             "loaded..."
#             sys.stdout.flush()
#
#             junctions_distrib = {}
#             junctions_list = {}
#
#             for (chrome,
#                  strand), int_seg_ids in segmentation._bisect_seg_ids.iteritems():
#                 int_sta_coords = segmentation._bisect_seg_sta_coords[
#                     (chrome, strand)]
#                 int_sto_coords = segmentation._bisect_seg_sto_coords[
#                     (chrome, strand)]
#
#                 prev_seg = segmentation.segments[int_seg_ids[0]]
#                 prev_sta = int_sta_coords[0]
#                 prev_sto = int_sto_coords[0]
#                 prev_len = prev_sto - prev_sta + 1
#
#                 for cur_seg_id, cur_sta, cur_sto in zip(int_seg_ids[1:],
#                                                         int_sta_coords[1:],
#                                                         int_sto_coords[1:]):
#                     cur_seg = segmentation.segments[cur_seg_id]
#                     assert cur_sta - prev_sto == 1
#                     cur_len = cur_sto - cur_sta + 1
#
#                     if strand == '+':
#                         up_seg_type = prev_seg.seg_type
#                         up_gene_name = prev_seg.gene_name
#                         up_gene_ID = prev_seg.gene_ID
#                         down_seg_type = cur_seg.seg_type
#                         down_gene_name = cur_seg.gene_name
#                         down_gene_ID = cur_seg.gene_ID
#                         jun_pos = cur_sta
#                         jun_up = min(max_flank, prev_len / 2 + prev_len % 2)
#                         jun_down = min(max_flank, cur_len - (
#                         cur_len / 2 + cur_len % 2) - 1)
#                     else:
#                         up_seg_type = cur_seg.seg_type
#                         up_gene_name = cur_seg.gene_name
#                         up_gene_ID = cur_seg.gene_ID
#                         down_seg_type = prev_seg.seg_type
#                         down_gene_name = prev_seg.gene_name
#                         down_gene_ID = prev_seg.gene_ID
#                         jun_pos = prev_sto
#                         jun_up = min(max_flank, cur_len / 2 + cur_len % 2)
#                         jun_down = min(max_flank, prev_len - (
#                         prev_len / 2 + prev_len % 2) - 1)
#
#                     # determine gene names associated with junction
#                     associated_gene_name = set()
#                     if up_gene_name <> 'inter' and up_gene_name <> 'telo':
#                         associated_gene_name.add(up_gene_name)
#                     if down_gene_name <> 'inter' and down_gene_name <> 'telo':
#                         associated_gene_name.add(down_gene_name)
#                     associated_gene_name = ", ".join(
#                         sorted(associated_gene_name))
#                     # and gene IDs
#                     associated_gene_ID = set()
#                     if up_gene_ID <> 'inter' and up_gene_ID <> 'telo':
#                         associated_gene_ID.add(up_gene_ID)
#                     if down_gene_ID <> 'inter' and down_gene_ID <> 'telo':
#                         associated_gene_ID.add(down_gene_ID)
#                     associated_gene_ID = ", ".join(sorted(associated_gene_ID))
#
#                     for mname, up_types, down_types, up_in_gene, down_in_gene, within_same_gene,
# minU, minD in maps_to_gen:
#                         if up_seg_type not in up_types or down_seg_type not in down_types: contin
# ue
#                         if within_same_gene and up_gene_ID <> down_gene_ID: continue
#                         if not (
#                         within_same_gene) and up_gene_ID == down_gene_ID: continue
#
#                         # check for min lengths of upstream and downstream
#                         if minU < -jun_up: continue
#                         if jun_down < minD: continue
#
#                         # positions jun_up and jun_down are included in the junction
#                         jun_rec = (chrome, strand, jun_pos, -jun_up, jun_down,
#                                    associated_gene_name, associated_gene_ID)
#                         junctions_list.setdefault(mname, []).append(jun_rec)
#                         # update junction frequency distribution for spanning region
#                         tmpd = junctions_distrib.setdefault(mname, {})
#                         for p in xrange(-jun_up, jun_down + 1):
#                             tmpd[p] = tmpd.get(p, 0) + 1
#
#                     prev_seg = cur_seg
#                     prev_sta = cur_sta
#                     prev_sto = cur_sto
#                     prev_len = cur_len
#
#             iCount.files.osutils.make_dir(
#                 _fname_RNAmap_folder(aname, annot_ver))
#             print
#             "saving data on %s maps" % len(junctions_list.keys())
#             sys.stdout.flush()
#             for mname, up_types, down_types, up_in_gene, down_in_gene, within_same_gene, minU, mi
# nD in maps_to_gen:
#                 jlist = junctions_list.get(mname, [])
#                 jlist.sort()
#                 all_junctions_cn = len(jlist)
#                 jdistrib = junctions_distrib.get(mname, {})
#
#                 maxU = min([0] + jdistrib.keys())
#                 maxD = max([0] + jdistrib.keys())
#                 junction_name = 'junction'
#                 desc_up = "regions of type: %s" % ", ".join(up_types)
#                 desc_down = "regions of type: %s" % ", ".join(down_types)
#
#                 fout = open(_fname_RNAmap(aname, annot_ver, mname), "wt")
#                 fout.write(
#                     '%s\t%s\t%s\n' % (mname, junction_name, all_junctions_cn))
#                 fout.write('%s\t%s\t%s\t%s\n' % (
#                 maxU, maxD, ['flanking', 'within'][up_in_gene],
#                 ['flanking', 'within'][down_in_gene]))
#                 fout.write('%s\n%s\n' % (desc_up, desc_down))
#                 for jun_rec in jlist:
#                     fout.write("%s\n" % "\t".join([str(x) for x in jun_rec]))
#                 fout.write('\n')
#
#                 vals = jdistrib.items()
#                 vals.sort()
#                 fout.write('%s\n' % "\t".join([str(p) for p, v in vals]))
#                 fout.write('%s\n' % "\t".join([str(v) for p, v in vals]))
#                 fout.close()
#             print
#
#
# ## BRANCH POINT data
# import sys
# import iCount
# import iCount.genomes
# import gzip
#
#
# def read_store(fname_in, assembly, annot_ver, map_name, junction_name, descU,
#                descD, fname_out, I_chrome, I_strand, I_start_pos, I_end_pos,
#                I_dist_to_intron_3ss, I_score):
#     maxU = -300
#     maxD = 100
#     minU = -100
#     minD = 20
#
#     fin = gzip.open(fname_in, "rt")
#     segmentation = iCount.genomes.annotation.Segmentation(assembly=assembly,
#                                                           annotation_version=annot_ver)
#     # read BP predictions and group them by containing intron
#     BP_by_intron = {}
#     all_introns = set()
#     for r in fin:
#         r = r.rstrip('\n\r').split('\t')
#         chrome = r[I_chrome]
#         intron_strand = r[I_strand]
#         intron_start_pos = int(r[I_start_pos]) - 1  # convert to zero-based
#         intron_end_pos = int(
#             r[I_end_pos]) - 1  # convert to zero-based, inclusive
#         assert (intron_start_pos <= intron_end_pos)
#         if r[I_dist_to_intron_3ss] == '.':
#             continue
#         dist_to_intron_3ss = int(r[I_dist_to_intron_3ss])
#         SVM_score = float(r[I_score])
#
#         assert (intron_start_pos <= intron_end_pos)
#         # take only the downstream half of intron
#         intron_half_len = intron_end_pos - intron_start_pos + 1
#         intron_half_len = intron_half_len / 2 + intron_half_len % 2
#         if intron_strand == '+':
#             BP_position = intron_end_pos - dist_to_intron_3ss
#             intron_start_pos = intron_end_pos - intron_half_len
#         else:
#             BP_position = intron_start_pos + dist_to_intron_3ss
#             intron_end_pos = intron_start_pos + intron_half_len
#         assert (intron_start_pos <= intron_end_pos)
#
#         if intron_start_pos <= BP_position <= intron_end_pos:
#             k = (chrome, intron_strand, intron_start_pos, intron_end_pos)
#             v = (SVM_score, BP_position)
#             all_introns.add(k)
#             if SVM_score > 0.0 and v not in BP_by_intron.setdefault(k, []):
#                 BP_by_intron.setdefault(k, []).append(v)
#
#     # sort and write
#     print
#     "##########################################################"
#     print
#     "fname_in:", fname_in
#     print
#     "fname_out:", fname_out
#     print
#     "all introns:", len(all_introns)
#     print
#     "introns with positive scores:", len(BP_by_intron)
#     print
#     # take best BP(s) from each intron, group by BP position. More than one BP can be choosen if
# equally well scored
#     print
#     "select top scored BP(s) in each intron"
#     by_pos = {}
#     for (chrome, intron_strand, intron_start_pos,
#          intron_end_pos), tmpl in BP_by_intron.iteritems():
#         tmpl.sort(reverse=True)  # by decreasing score
#         top_score = tmpl[0][0]
#         tmpl = [BP_position for SVM_score, BP_position in tmpl if
#                 SVM_score >= top_score]
#         if len(tmpl) > 1:
#             print
#             "multiple best score BPs:", chrome, intron_strand, intron_start_pos, intron_end_pos
#             if intron_strand == '+':
#                 tmpl = [(intron_end_pos - BP_position, BP_position) for
#                         BP_position in tmpl]
#             else:
#                 tmpl = [(BP_position - intron_start_pos, BP_position) for
#                         BP_position in tmpl]
#             tmpl.sort()  # by increasing downstream from BP to intron end
#             for x in tmpl:
#                 print
#                 x
#             BL_position = tmpl[0][1]
#             print
#         else:
#             BP_position = tmpl[0]
#
#         v = (intron_start_pos, intron_end_pos)
#         tmpl = by_pos.setdefault((chrome, intron_strand, BP_position), [])
#         if v not in tmpl:
#             tmpl.append(v)
#     print
#
#     # in case more than one intron spans BP, select largest intron
#     junctions_distrib = {}
#     print
#     "checking for multiple introns spanning same BP"
#     junctions_list = []
#     for (chrome, intron_strand, BP_position), tmpl in by_pos.iteritems():
#         if len(tmpl) > 1:
#             print
#             "multiple introns for BP:", chrome, intron_strand, BP_position
#             # select span that covers all
#             for x in tmpl:
#                 print
#                 x
#             tmpl_s, tmpl_e = zip(*tmpl)
#             intron_start_pos = min(tmpl_s)
#             intron_end_pos = max(tmpl_e)
#             #            tmpl = [(abs(intron_start_pos-intron_end_pos), intron_start_pos, intron_
# end_pos) for intron_start_pos, intron_end_pos in tmpl]
#             #            tmpl.sort(reverse=True)
#             #            intron_start_pos, intron_end_pos = tmpl[0][1:]
#             print
#             intron_start_pos, intron_end_pos
#             print
#         else:
#             intron_start_pos, intron_end_pos = tmpl[0]
#
#         if intron_strand == '+':
#             spanU = intron_start_pos - BP_position
#             spanD = intron_end_pos - BP_position
#         else:
#             spanU = BP_position - intron_end_pos
#             spanD = BP_position - intron_start_pos
#         assert (spanU <= 0 and spanD >= 0)
#
#         if minU < spanU: continue
#         if spanD < minD: continue
#
#         spanU = max(maxU, spanU)
#         spanD = min(maxD, spanD)
#         for p in xrange(spanU, spanD + 1):
#             junctions_distrib[p] = junctions_distrib.get(p, 0) + 1
#
#         _, _, hit_int_same, hit_int_anti = segmentation.get_annotation(chrome,
#                                                                        intron_strand,
#                                                                        BP_position)
#         if hit_int_anti:
#             associated_gene_name_anti = hit_int_anti.segment.gene_name
#         else:
#             associated_gene_name_anti = '?'
#
#         if hit_int_same:
#             associated_gene_name_same = hit_int_same.segment.gene_name
#             if (
#                     associated_gene_name_same == 'inter' or associated_gene_name_same == 'telo')
# and associated_gene_name_anti <> '?':
#                 associated_gene_name = associated_gene_name_anti
#             else:
#                 associated_gene_name = associated_gene_name_same
#         else:
#             associated_gene_name = associated_gene_name_anti
#
#         junctions_list.append((
#                               chrome, intron_strand, BP_position, spanU, spanD,
#                               associated_gene_name))
#
#     junctions_list.sort()
#
#     # adjust maxU and maxD to actual span (if shorter than originally given)
#     maxU = max(maxU, min([0] + junctions_distrib.keys()))
#     maxD = min(maxD, max([0] + junctions_distrib.keys()))
#     fout = open(fname_out, "wt")
#     all_junctions_cn = len(junctions_list)
#     fout.write('%s\t%s\t%s\n' % (map_name, junction_name, all_junctions_cn))
#     fout.write('%s\t%s\t%s\t%s\n' % (maxU, maxD, 'within', 'within'))
#     fout.write('%s\n%s\n' % (descU, descD))
#     for r in junctions_list:
#         ostr = "\t".join([str(x) for x in r])
#         fout.write("%s\n" % ostr)
#     fout.write('\n')
#
#     vals = junctions_distrib.items()
#     vals.sort()
#     fout.write('%s\n' % "\t".join([str(p) for p, v in vals]))
#     fout.write('%s\n' % "\t".join([str(v) for p, v in vals]))
#     fout.close()
#
#
# """
# ### *********************
# ### mm9 branch points
# fname_in = "prep/Mmus.BPpred.gz"
# assembly = "mm9"
# annot_ver = "ensembl59"
# map_name = "BP best in intron"
# junction_name = "BP"
# descU = 'usptream of branch point'
# descD = 'downstream of branch point'
# fname_out = "%s_%s_RNAmap_BP_best_one_in_intron_v0.tab" % (assembly, annot_ver)
# read_store(fname_in, assembly, annot_ver, map_name, junction_name, descU, descD, fname_out, 1, 2,
#  3, 4, 6, 13)
#
# ### *********************
# ### hg18 branch points
# fname_in = "prep/hg19_remapped_Hsap.BP.predictions.gz"
# assembly = "hg19"
# annot_ver = "ensembl59"
# map_name = "BP best in intron"
# junction_name = "BP"
# descU = 'usptream of branch point'
# descD = 'downstream of branch point'
# fname_out = "%s_%s_RNAmap_BP_best_one_in_intron_v0.tab" % (assembly, annot_ver)
# read_store(fname_in, assembly, annot_ver, map_name, junction_name, descU, descD, fname_out, 2, 3,
#  4, 5, 8, 16)
#
# """
