"""
k-mer enrichment analysis
-------------------------

Read bedGraph with cross-linked sites.
Count k-mer frequencies.
Perform permutation analysis to determine significance of observed k-mer
frequencies.
Return ranked list of k-mer enrichment.

"""



# description and parameters needed for the analysis
analysis_name = 'kmers'
analysis_description_short = 'k-mer enrichment analysis'
analysis_description = 'Determine enriched k-mers in vicinity of ' \
                       'cross-linked sites.'

params_opt = [
    (
        'k', 'int_range', (5, 2, 9), False,
        'Length of k-mers.'
    ),
    (
        'rnd', 'int_range', (100, 10, 250), False,
        'Number of random permutations needed to determine statistical '
        'significance.'
    ),
    (
        'intervals', 'intervals', ([(-30, -10), (10, 30)], -50, 50), False,
        'Count k-mers in intervals relative to cross-link.',
    ),
    (
        'report_region', 'intervals', ([(-50, 50)], -100, 100), False,
        'Report frequency of k-mer occurence at each position in region.'
    ),
    (
        'regions', 'str_list', [], False,
        'Consider only cross-links in regions of selected types.'
    ),
    (
        'chromosomes', 'str_list', [], False,
        'Consider only cross-links on listed chromosomes.'

    ),
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
        'enrichment', 'tsv', 'out',
        '(output) Tab-delimited file with information on enriched k-mers.'
    ),
]


# def gen_mer(s):
#     nucs = ['A', 'T', 'C', 'G']
#     if s == 0:
#         return []
#     if s == 1:
#         return [x for x in nucs]
#     sub_mers = gen_mer(s-1)
#
#     ret_mers = []
#     for m in sub_mers:
#         for n in nucs:
#             ret_mers.append( m+n)
#     return ret_mers
#
# def score_kmers(bed_d, segmentation, random_perms, k, intervals, (report_for_region_start, report_for_region_stop), valid_region_types, selected_chromes, analysis_id=None):
#     """ to preserve memory, bed_d will be destroyed in the process """
#     random.seed(42)
#
#     all_Kmers = gen_mer(k)
#     all_Kmers_expected = 4**k
#     assert(len(set(all_Kmers)) == all_Kmers_expected)
#
#     undef_z_score = 100000.0
#     tot_lens = [interval[1] - interval[0] + 1 for interval in intervals]
#     interval_spans = [tot_len - (k-1) for tot_len in tot_lens]
#     interval_str = ".".join(["%s..%s" % interval for interval in intervals])
#
# # k = 5, k_c = 5/2 = 2
# # kkkkk                    kkkkk (margins)
# #   01234567890123456789012345   (user defined region)
# #
# # k = 2, k_c = 4/2 - 1 = 1
# # on +:
# #  kkkk                     kkkk
# #   01234567890123456789012345   (user defined region)
# #
# # on -:
# # kkkk                     kkkk
# #   54321098765432109876543210   (user defined region)
# #
#
# #0123456789
# #       ---
#     kmer_center_pos = k / 2 # center of kmer if odd length
#     if k % 2 == 0:
#         kmer_center_pos -= 1 # one position prior to center if even length
#         report_for_region_start_with_margin = report_for_region_start - kmer_center_pos
#         report_for_region_stop_with_margin = report_for_region_stop + kmer_center_pos + 1
#     else:
#         report_for_region_start_with_margin = report_for_region_start - kmer_center_pos
#         report_for_region_stop_with_margin = report_for_region_stop + kmer_center_pos
#
#     tot_len_rep = report_for_region_stop - report_for_region_start + 1
#     tot_len_rep_with_margin = report_for_region_stop_with_margin - report_for_region_start_with_margin + 1
#     interval_span_rep = tot_len_rep_with_margin - (k-1)
#
#     strands = ['-', '+']
#     xlinks_kmer_counts = {}
#     xlinks_kmer_counts_by_pos = {} # iBind
#     xlinks_kmer_counts_rnd = {}
#     count_kmer_counts = {}
#     count_kmer_counts_by_pos = {} # iBind
#     count_kmer_counts_rnd = {}
#     xlinks_kmer_occurrences = {}
#     xlinks_kmer_occurrences_rnd = {}
#     count_kmer_occurrences = {}
#     count_kmer_occurrences_rnd = {}
#     for kmer in all_Kmers:
#         count_kmer_counts_rnd[kmer] = [0.0]*random_perms
#         xlinks_kmer_counts_rnd[kmer] = [0.0]*random_perms
#         count_kmer_occurrences_rnd[kmer] = [0.0]*random_perms
#         xlinks_kmer_occurrences_rnd[kmer] = [0.0]*random_perms
#         xlinks_kmer_counts_by_pos[kmer] = [0.0]*interval_span_rep
#         count_kmer_counts_by_pos[kmer] = [0.0]*interval_span_rep
#
#     # list of gene_name, position, frequency of kmers in vicinity
#     bygene_bypos_bykmer = {}
#
#     xlinks_total = {}
#     xlinks_filtered_out = {}
#     count_total = {}
#     count_filtered_out = {}
#     prev_chrome = None
#     prev_prog = 5
#     tot = len(bed_d)
#     while bed_d:
#         (chrome, strand), tmpd = bed_d.popitem()
#         if selected_chromes <> 'whole genome' and chrome not in selected_chromes:
#             continue
#
#         if len(tmpd) == 0:
#             print "Warning, no hits on:", chrome, strand
#             continue
#         if chrome <> prev_chrome:
#             chrome_seq = iCount.genomes.assembly[segmentation.assembly].get_chromosome(chrome)
#             valid_nucs = set(['A', 'T', 'G', 'C'])
#             chrome_seq = chrome_seq.upper()
#             chrome_seq = "".join([n if n in valid_nucs else 'N' for n in chrome_seq])
#             prev_chrome = chrome
#
#         for (p_start, p_end), value in tmpd.iteritems():
#             # (ii) iCLIP reads antisense to the transcriptional direction of the associated gene, and reads that mapped to non-annotated genomic regions, were removed before proceeding to further analysis.
#             hit_int, orient, _, _ = segmentation.get_annotation(chrome, strand, p_start)
#             if not hit_int:
#                 xlinks_total['undefined'] = xlinks_total.get('undefined', 0) + 1
#                 xlinks_filtered_out['undefined'] = xlinks_filtered_out.get('undefined', 0) + 1
#                 count_total['undefined'] = count_total.get('undefined', 0) + value
#                 count_filtered_out['undefined'] = count_filtered_out.get('undefined', 0) + value
#                 continue
#
#             seg = hit_int.segment
#             gene_name = seg.gene_name
#             rtype = seg.seg_type
#             rtype_biotype = seg.biotype
#             rtype_k = "%s, %s" % (rtype, rtype_biotype)
#             xlinks_total[rtype_k] = xlinks_total.get(rtype_k, 0) + 1
#             count_total[rtype_k] = count_total.get(rtype_k, 0) + value
#
#             if rtype == 'inter' or rtype == 'telo':
#                 gene_name = ''
#
#             if orient <> 'same':
#                 xlinks_filtered_out[rtype_k] = xlinks_filtered_out.get(rtype_k, 0) + 1
#                 count_filtered_out[rtype_k] = count_filtered_out.get(rtype_k, 0) + value
#                 continue
#
#             for valid_rtype in valid_region_types:
#                 if type(valid_rtype) == str:
#                     if valid_rtype == rtype:
#                         break
#                 elif type(valid_rtype) == tuple:
#                     if valid_rtype == (rtype, rtype_biotype):
#                         break
#                 else:
#                     assert(False)
#             else:
#                 xlinks_filtered_out[rtype_k] = xlinks_filtered_out.get(rtype_k, 0) + 1
#                 count_filtered_out[rtype_k] = count_filtered_out.get(rtype_k, 0) + value
#                 continue
#
#             # count occurence of kmers in report_for_region interval (for "iBind" matrix)
#             if strand == '+':
#                 sta_pos = p_start + report_for_region_start_with_margin
#                 sto_pos = p_start + report_for_region_stop_with_margin
#             else:
#                 sta_pos = p_start - report_for_region_stop_with_margin
#                 sto_pos = p_start - report_for_region_start_with_margin
#             genome_seq = chrome_seq[sta_pos:sto_pos+1]
#             if len(genome_seq) != tot_len_rep_with_margin:
#                 print "problem with filling kmer positional occurence matrix:", chrome, p_start, value, sta_pos, sto_pos
#                 print "it is %s, should be %s" % (len(genome_seq), tot_len_rep_with_margin)
#                 print genome_seq
#                 print ""
#             else:
#                 if strand == '-':
#                     genome_seq = iCount.misclibs.reverse_complement(genome_seq)
#
#                 for i in range(interval_span_rep):
#                     kmer = genome_seq[i:i+k]
#                     if 'N' in kmer: continue
#                     xlinks_kmer_counts_by_pos[kmer][i] += 1
#                     count_kmer_counts_by_pos[kmer][i] += value
#
#             if gene_name:
#                 tmpd_bygene = bygene_bypos_bykmer.setdefault(gene_name, {})
#                 tmpd_bypos = tmpd_bygene.setdefault((chrome, strand, p_start, value, seg), {})
#             else:
#                 tmpd_bypos = {} # these x-links will be ignored because they do not belong to a gene
#
#             # count kmers for significance testing
#             for interval, tot_len, interval_span in zip(intervals, tot_lens, interval_spans):
#                 # extract genome sequence
#                 if strand == '+':
#                     sta_pos = p_start + interval[0]
#                     sto_pos = p_start + interval[1]
#                 else:
#                     sta_pos = p_start - interval[1]
#                     sto_pos = p_start - interval[0]
#
#                 genome_seq = chrome_seq[sta_pos:sto_pos+1]
#                 if len(genome_seq) != tot_len:
#                     print "problems with:", chrome, p_start, value
#                     print "              ", sta_pos, sto_pos
#                     print genome_seq
#                     continue
#
#                 if strand == '-':
#                     genome_seq = iCount.misclibs.reverse_complement(genome_seq)
#
#                 # calculate on true data
#                 ps = [genome_seq[i:i+k] for i in range(interval_span)]
#                 for kmer in ps:
#                     count_kmer_counts[kmer] = count_kmer_counts.get(kmer, 0) + value
#                     xlinks_kmer_counts[kmer] = xlinks_kmer_counts.get(kmer, 0) + 1
#                     tmpd_bypos[kmer] = tmpd_bypos.get(kmer, 0) + 1
#                 for kmer in set(ps): # count each pentamer only once!
#                     count_kmer_occurrences[kmer] = count_kmer_occurrences.get(kmer, 0) + value
#                     xlinks_kmer_occurrences[kmer] = xlinks_kmer_occurrences.get(kmer, 0) + 1
#
#                 # randomly permute
#                 for r in range(random_perms):
#                     # pick random position within current segment
#                     rnd_pos = seg.random_position()
#
#                     if strand == '+':
#                         sta_pos = rnd_pos + interval[0]
#                         sto_pos = rnd_pos + interval[1]
#                     else:
#                         sta_pos = rnd_pos - interval[1]
#                         sto_pos = rnd_pos - interval[0]
#                     genome_seq = chrome_seq[sta_pos:sto_pos+1]
#                     if len(genome_seq) != tot_len:
#                         print "problem with:", chrome, p_start, value
#                         print "             ", sta_pos, sto_pos
#                         print genome_seq
#                         continue
#
#                     if strand == '-':
#                         genome_seq = iCount.misclibs.reverse_complement(genome_seq)
#
#                     ps = [genome_seq[i:i+k] for i in range(interval_span)]
#                     for kmer in ps:
#                         if 'N' in kmer: continue
#                         count_kmer_counts_rnd[kmer][r] += value
#                         xlinks_kmer_counts_rnd[kmer][r] += 1.0
#                     for kmer in set(ps): #count each pentamer only once!
#                         if 'N' in kmer: continue
#                         count_kmer_occurrences_rnd[kmer][r] += value
#                         xlinks_kmer_occurrences_rnd[kmer][r] += 1.0
#         if analysis_id is not None:
#             prog = 85 - (100000*len(bed_d)/tot)/1250
#             if prog > prev_prog:
#                 db.analysis_set_status(analysis_id, 'processing %s%%' % prog)
#                 prev_prog = prog
#
#     all_zero_by_pos = [0.0]*interval_span_rep
#     res = {}
#     tot = len(all_Kmers)
#     for progi, kmer in enumerate(all_Kmers):
#         # count (count occurences)
#         count_v = count_kmer_counts.get(kmer, 0.0)
#         count_cns = count_kmer_counts_rnd[kmer]
#         assert( len(count_cns) == random_perms)
#         count_st = numpy.std(count_cns)
#         count_m = numpy.mean(count_cns)
#         try:
#             count_z = (count_v - count_m) / count_st
#         except:
#             count_z = undef_z_score
#         count_v_pval_ge = float(sum([int(cn >= count_v) for cn in count_cns]))/max(1.0, len(count_cns))
#         count_v_pval_le = float(sum([int(cn <= count_v) for cn in count_cns]))/max(1.0, len(count_cns))
#         count_v_pval = min(count_v_pval_ge, count_v_pval_le)
#
#         # x-links (count occurences)
#         xlinks_v = xlinks_kmer_counts.get(kmer, 0.0)
#         xlinks_cns = xlinks_kmer_counts_rnd[kmer]
#         assert( len(xlinks_cns) == random_perms)
#         xlinks_st = numpy.std(xlinks_cns)
#         xlinks_m = numpy.mean(xlinks_cns)
#         try:
#             xlinks_z = (xlinks_v - xlinks_m) / xlinks_st
#         except:
#             xlinks_z = undef_z_score
#         xlinks_v_pval_ge = float(sum([int(cn >= xlinks_v) for cn in xlinks_cns]))/max(1.0, len(xlinks_cns))
#         xlinks_v_pval_le = float(sum([int(cn <= xlinks_v) for cn in xlinks_cns]))/max(1.0, len(xlinks_cns))
#         xlinks_v_pval = min(xlinks_v_pval_ge, xlinks_v_pval_le)
#
#         # count (single count for each kmer occurrence)
#         occur_count_v = count_kmer_occurrences.get(kmer, 0.0)
#         occur_count_cns = count_kmer_occurrences_rnd[kmer]
#         assert( len(occur_count_cns) == random_perms)
#         occur_count_st = numpy.std(occur_count_cns)
#         occur_count_m = numpy.mean(occur_count_cns)
#         try:
#             occur_count_z = (occur_count_v - occur_count_m) / occur_count_st
#         except:
#             occur_count_z = undef_z_score
#         occur_count_v_pval_ge = float(sum([int(cn >= occur_count_v) for cn in occur_count_cns]))/max(1.0, len(occur_count_cns))
#         occur_count_v_pval_le = float(sum([int(cn <= occur_count_v) for cn in occur_count_cns]))/max(1.0, len(occur_count_cns))
#         occur_count_v_pval = min(occur_count_v_pval_ge, occur_count_v_pval_le)
#
#         # x-links (single count for each kmer occurence)
#         occur_xlinks_v = xlinks_kmer_occurrences.get(kmer, 0.0)
#         occur_xlinks_cns = xlinks_kmer_occurrences_rnd[kmer]
#         assert( len(occur_xlinks_cns) == random_perms)
#         occur_xlinks_st = numpy.std(occur_xlinks_cns)
#         occur_xlinks_m = numpy.mean(occur_xlinks_cns)
#         try:
#             occur_xlinks_z = (occur_xlinks_v - occur_xlinks_m) / occur_xlinks_st
#         except:
#             occur_xlinks_z = undef_z_score
#         occur_xlinks_v_pval_ge = float(sum([int(cn >= occur_xlinks_v) for cn in occur_xlinks_cns]))/max(1.0, len(occur_xlinks_cns))
#         occur_xlinks_v_pval_le = float(sum([int(cn <= occur_xlinks_v) for cn in occur_xlinks_cns]))/max(1.0, len(occur_xlinks_cns))
#         occur_xlinks_v_pval = min(occur_xlinks_v_pval_ge, occur_xlinks_v_pval_le)
#
#         # normalize _by_pos statistics
#         # count_kmer_counts_by_pos
#         if kmer in count_kmer_counts_by_pos:
#             tmpl = count_kmer_counts_by_pos[kmer]
#             norm_count_kmer_counts_by_pos = [v/count_m for v in tmpl]
# #            norm_count_kmer_counts_by_pos = [v/1.0 for v in tmpl]
#         else:
#             norm_count_kmer_counts_by_pos = all_zero_by_pos
#         # normalize _by_pos statistics
#         # xlinks_kmer_counts_by_pos
#         if kmer in xlinks_kmer_counts_by_pos:
#             tmpl = xlinks_kmer_counts_by_pos[kmer]
#             norm_xlinks_kmer_counts_by_pos = [v/xlinks_m for v in tmpl]
# #            norm_xlinks_kmer_counts_by_pos = [v/1.0 for v in tmpl]
#         else:
#             norm_xlinks_kmer_counts_by_pos = all_zero_by_pos
#
#         res[kmer] = (
#             count_z, count_v, count_v_pval, count_m, count_st, count_cns,
#             xlinks_z, xlinks_v, xlinks_v_pval, xlinks_m, xlinks_st, xlinks_cns,
#             occur_count_z, occur_count_v, occur_count_v_pval, occur_count_m, occur_count_st, occur_count_cns,
#             occur_xlinks_z, occur_xlinks_v, occur_xlinks_v_pval, occur_xlinks_m, occur_xlinks_st, occur_xlinks_cns,
#             norm_count_kmer_counts_by_pos, norm_xlinks_kmer_counts_by_pos
#         )
#
#         if analysis_id is not None:
#             prog = 85 + (100000*(progi+1)/tot)/10000
#             if prog > prev_prog:
#                 db.analysis_set_status(analysis_id, 'processing %s%%' % prog)
#                 prev_prog = prog
#     return (res, kmer_center_pos), xlinks_total, xlinks_filtered_out, count_total, count_filtered_out, bygene_bypos_bykmer
#
# class KmerScores:
#     def __init__(self, kmer, xlinks_z, xlinks_v, xlinks_pval, xlinks_rnd_mean, xlinks_rnd_stdev, xlinks_rnd_vs, count_z, count_v, count_pval, count_rnd_mean, count_rnd_stdev, count_rnd_vs):
#         self.kmer = kmer
#
#         self.xlinks_z = xlinks_z
#         self.xlinks_v = xlinks_v
#         self.xlinks_pval = xlinks_pval
#         self.xlinks_rnd_mean = xlinks_rnd_mean
#         self.xlinks_rnd_stdev = xlinks_rnd_stdev
#         self.xlinks_rnd_vs = xlinks_rnd_vs
#
#         self.count_z = count_z
#         self.count_v = count_v
#         self.count_pval = count_pval
#         self.count_rnd_mean = count_rnd_mean
#         self.count_rnd_stdev = count_rnd_stdev
#         self.count_rnd_vs = count_rnd_vs
#
# def read_results(in_fname):
#     res = {}
#     f = iCount.nc_open(in_fname, "rt")
#     r = f.readline()
#     while r.startswith('#'):
#         r = f.readline()
#     assert(r.startswith('kmer'))
#     r = f.readline()
#     while r:
#         r = r.rstrip('\n\r').split('\t')
#         kmer, xlinks_z, count_z, xlinks_v, xlinks_pval, count_v, count_pval, xlinks_rnd_mean, count_rnd_mean, xlinks_rnd_stdev, count_rnd_stdev, xlinks_rnd_vs, count_rnd_vs = r
#         xlinks_z = float(xlinks_z)
#         count_z = float(count_z)
#         xlinks_v = float(xlinks_v)
#         xlinks_pval = float(xlinks_pval)
#         count_v = float(count_v)
#         count_pval = float(count_pval)
#         xlinks_rnd_mean = float(xlinks_rnd_mean)
#         count_rnd_mean = float(count_rnd_mean)
#         xlinks_rnd_stdev = float(xlinks_rnd_stdev)
#         count_rnd_stdev = float(count_rnd_stdev)
#         xlinks_rnd_vs = [float(x) for x in xlinks_rnd_vs.split(',')]
#         count_rnd_vs = [float(x) for x in count_rnd_vs.split(',')]
#
#         res[kmer] = KmerScores(kmer, xlinks_z, xlinks_v, xlinks_pval, xlinks_rnd_mean, xlinks_rnd_stdev, xlinks_rnd_vs, count_z, count_v, count_pval, count_rnd_mean, count_rnd_stdev, count_rnd_vs)
#
#         r = f.readline()
#     return res
#
# #######
# # determine filename from parameters
# def process_file(in_fname, mapped_to, annot_ver, out_fname_pref, random_perms, k, intervals, report_for_region, valid_region_types, selected_chromes, segmentation=None, analysis_id=None):
#     print "*********************************"
#     print time.ctime()
#     print "scoring kmers in file:", in_fname
#     print "saving to:", out_fname_pref
#     print "parameters:"
#     print "\trandom permutations:", random_perms
#     print "\tk:", k
#     print "\tintervals:", intervals
#     print "\treport_for_region:", str(report_for_region)
#     print "\tvalid region types:", valid_region_types
#     print "\tselected chromes:", selected_chromes
#     sys.stdout.flush()
#     if analysis_id is not None:
#         db.analysis_set_status(analysis_id, 'processing started')
#         prev_prog = 0
#
#     bed_d = BED.load_bedGraph(in_fname)
#     bed_d, bed_h = BED.load_bedGraph(in_fname, return_tracks_descriptions=True)
#     assert(len(bed_h) == 1)
#     bed_h = bed_h[0]
#     try:
#         bed_res_type = bed_h.split('res_type=')[1].split(' ')[0].strip('"')
#     except:
#         bed_res_type = 'unknown_counts'
#     bed_res_type_desc = beds.desc_by_result_type.get(bed_res_type, None)
#     if bed_res_type_desc is not None:
#         bed_res_type = bed_res_type_desc[-1]
#     print "input res_type: %s" % bed_res_type
#     sys.stdout.flush()
#
#     if segmentation is None or segmentation.assembly <> mapped_to or segmentation.annotation_version <> annot_ver:
#         segmentation = iCount.genomes.annotation.Segmentation(assembly=mapped_to, annotation_version=annot_ver)
#         print "segmentation loaded"
#         sys.stdout.flush()
#
#     if analysis_id is not None:
#         db.analysis_set_status(analysis_id, 'processing 5%')
#         prev_prog = 5
#
#     (stats, kmer_center_pos), xlinks_total, xlinks_filtered_out, count_total, count_filtered_out, bygene_bypos_bykmer = score_kmers(bed_d, segmentation, random_perms, k, intervals, report_for_region, valid_region_types, selected_chromes, analysis_id)
#     print "kmer statistics calculated"
#     sys.stdout.flush()
#
#     # write stats - count all occurrences of a motif within the defined interval(s) relative to each x-link
#     fout = iCount.nc_open(out_fname_pref+"_countAll.tab.gz", "wt")
#     o_stats = []
#     for (kmer, (count_z,  count_v,  count_v_pval,  count_m,  count_st,  count_cns,
#                 xlinks_z, xlinks_v, xlinks_v_pval, xlinks_m, xlinks_st, xlinks_cns,
#                 occur_count_z,  occur_count_v,  occur_count_v_pval,  occur_count_m,  occur_count_st,  occur_count_cns,
#                 occur_xlinks_z, occur_xlinks_v, occur_xlinks_v_pval, occur_xlinks_m, occur_xlinks_st, occur_xlinks_cns,
#                 norm_count_kmer_counts_by_pos, norm_xlinks_kmer_counts_by_pos)) in stats.iteritems():
#         o_stats.append((xlinks_z + count_z,
#                        xlinks_z, count_z, occur_xlinks_z, occur_count_z,
#                        xlinks_v, xlinks_v_pval,
#                        count_v, count_v_pval,
#                        occur_xlinks_v, occur_xlinks_v_pval,
#                        occur_count_v, occur_count_v_pval,
#                        xlinks_m, count_m, occur_xlinks_m, occur_count_m,
#                        xlinks_st, count_st, occur_xlinks_st, occur_count_st,
#                        xlinks_cns, count_cns, occur_xlinks_cns, occur_count_cns,
#                        norm_count_kmer_counts_by_pos, norm_xlinks_kmer_counts_by_pos,
#                        kmer))
#     o_stats.sort(reverse=True)
#     fout.write("# All occurrences of a kmer in interval(s) (%s) relative to each x-link are counted.\n" % intervals)
#     fout.write("# x-link sites - each occurrence is weighted by 1.0\n")
#     fout.write("# %s - each occurrence is weighted by this score at the x-link.\n" % bed_res_type)
#     fout.write("# Random reference data was generated %s times by random shuffling of iCLIP x-link positions within corresponding genome segments (within same genes).\n" % random_perms)
#     fout.write("kmer\tz-score [x-link sites]\tz-score [%s]\ttrue score [x-link sites]\ttrue score p.val [x-link sites]\ttrue score [%s]\ttrue score p.val [%s]\tmean random score [x-link sites]\tmean random score [%s]\tstdev random score [x-link sites]\tstdev random score [%s]\trandom scores [x-link sites]\trandom scores [%s]\n" % ((bed_res_type,)*6))
#     for (zs, xlinks_z, count_z, occur_xlinks_z, occur_count_z,
#                        xlinks_v, xlinks_v_pval,
#                        count_v, count_v_pval,
#                        occur_xlinks_v, occur_xlinks_v_pval,
#                        occur_count_v, occur_count_v_pval,
#                        xlinks_m, count_m, occur_xlinks_m, occur_count_m,
#                        xlinks_st, count_st, occur_xlinks_st, occur_count_st,
#                        xlinks_cns, count_cns, occur_xlinks_cns, occur_count_cns,
#                        norm_count_kmer_counts_by_pos, norm_xlinks_kmer_counts_by_pos,
#                        kmer) in o_stats:
# #        xlinks_cns = count_cns = ['intentionally blank']
#         fout.write("%s\t%0.3f\t%0.3f\t%0.1f\t%0.3g\t%0.1f\t%0.3g\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%s\t%s\n" % (kmer, xlinks_z, count_z, xlinks_v, xlinks_v_pval, count_v, count_v_pval, xlinks_m, count_m, xlinks_st, count_st, ",".join([str(x) for x in xlinks_cns]), ",".join([str(x) for x in count_cns])))
#     fout.write("\n")
#
#     # positional distribution
#     report_for_region_start, report_for_region_stop = report_for_region
#     fout.write("Positional distribution of k-mers in region %s..%s relative to x-link, normalized by mean random [x-link sites] score.\n" % (report_for_region_start, report_for_region_stop))
#     fout.write("Normalized frequencies of kmer centered at reported position(s) (i.e., position %s in kmer, zero-based indexing).\n" % (kmer_center_pos))
#     fout.write("kmer\tz-score [x-link sites]\tmean random score [x-link sites]")
#     for p in range(report_for_region_start, report_for_region_stop+1):
#         fout.write("\t%s" % p)
#     fout.write("\n")
#     for (zs, xlinks_z, count_z, occur_xlinks_z, occur_count_z,
#                        xlinks_v, xlinks_v_pval,
#                        count_v, count_v_pval,
#                        occur_xlinks_v, occur_xlinks_v_pval,
#                        occur_count_v, occur_count_v_pval,
#                        xlinks_m, count_m, occur_xlinks_m, occur_count_m,
#                        xlinks_st, count_st, occur_xlinks_st, occur_count_st,
#                        xlinks_cns, count_cns, occur_xlinks_cns, occur_count_cns,
#                        norm_count_kmer_counts_by_pos, norm_xlinks_kmer_counts_by_pos,
#                        kmer) in o_stats:
# #        xlinks_cns = count_cns = ['intentionally blank']
#         ostr = "\t".join(["%g" % v for v in norm_xlinks_kmer_counts_by_pos])
#         fout.write("%s\t%0.3f\t%0.3f\t%s\n" % (kmer, xlinks_z, xlinks_m, ostr))
#     fout.write("\n")
#
#     fout.write("Positional distribution of k-mers in region %s..%s relative to x-link, normalized by mean random [%s] score.\n" % (report_for_region_start, report_for_region_stop, bed_res_type))
#     fout.write("Normalized frequencies of kmer centered at reported position(s) (i.e., position %s in kmer, zero-based indexing).\n" % (kmer_center_pos))
#     fout.write("kmer\tz-score [%s]\tmean random score [%s]" % (bed_res_type, bed_res_type))
#     for p in range(report_for_region_start, report_for_region_stop+1):
#         fout.write("\t%s" % p)
#     fout.write("\n")
#     for (zs, xlinks_z, count_z, occur_xlinks_z, occur_count_z,
#                        xlinks_v, xlinks_v_pval,
#                        count_v, count_v_pval,
#                        occur_xlinks_v, occur_xlinks_v_pval,
#                        occur_count_v, occur_count_v_pval,
#                        xlinks_m, count_m, occur_xlinks_m, occur_count_m,
#                        xlinks_st, count_st, occur_xlinks_st, occur_count_st,
#                        xlinks_cns, count_cns, occur_xlinks_cns, occur_count_cns,
#                        norm_count_kmer_counts_by_pos, norm_xlinks_kmer_counts_by_pos,
#                        kmer) in o_stats:
# #        xlinks_cns = count_cns = ['intentionally blank']
#         ostr = "\t".join(["%g" % v for v in norm_count_kmer_counts_by_pos])
#         fout.write("%s\t%0.3f\t%0.3f\t%s\n" % (kmer, count_z, count_m, ostr))
#     fout.write("\n")
#     fout.close()
#     print "stats for all occurrences stored to:", out_fname_pref+"_countAll.tab.gz"
#     sys.stdout.flush()
#
#     # write stats - count only one (first) occurrence of a motif within the defined interval(s) relative to each x-link
#     fout = iCount.nc_open(out_fname_pref+"_countOne.tab.gz", "wt")
#     o_stats = []
#     for (kmer, (count_z,  count_v,  count_v_pval,  count_m,  count_st,  count_cns,
#                 xlinks_z, xlinks_v, xlinks_v_pval, xlinks_m, xlinks_st, xlinks_cns,
#                 occur_count_z,  occur_count_v,  occur_count_v_pval,  occur_count_m,  occur_count_st,  occur_count_cns,
#                 occur_xlinks_z, occur_xlinks_v, occur_xlinks_v_pval, occur_xlinks_m, occur_xlinks_st, occur_xlinks_cns,
#                 norm_count_kmer_counts_by_pos, norm_xlinks_kmer_counts_by_pos)) in stats.iteritems():
#         o_stats.append((occur_count_z + occur_xlinks_z,
#                        xlinks_z, count_z, occur_xlinks_z, occur_count_z,
#                        xlinks_v, xlinks_v_pval,
#                        count_v, count_v_pval,
#                        occur_xlinks_v, occur_xlinks_v_pval,
#                        occur_count_v, occur_count_v_pval,
#                        xlinks_m, count_m, occur_xlinks_m, occur_count_m,
#                        xlinks_st, count_st, occur_xlinks_st, occur_count_st,
#                        xlinks_cns, count_cns, occur_xlinks_cns, occur_count_cns,
#                        kmer))
#     o_stats.sort(reverse=True)
#     fout.write("# Only one occurrence of a kmer in interval(s) (%s) relative to each x-link is counted.\n" % intervals)
#     fout.write("# x-link sites - each occurrence is weighted by 1.0\n")
#     fout.write("# %s - each occurrence is weighted by this score at the x-link.\n" % bed_res_type)
#     fout.write("# Random reference data was generated %s times by random shuffling of iCLIP x-link positions within corresponding genome segments (within same genes).\n" % random_perms)
#     fout.write("kmer\tz-score [x-link sites]\tz-score [%s]\ttrue score [x-link sites]\ttrue score p.val [x-link sites]\ttrue score [%s]\ttrue score p.val [%s]\tmean random score [x-link sites]\tmean random score [%s]\tstdev random score [x-link sites]\tstdev random score [%s]\trandom scores [x-link sites]\trandom scores [%s]\n" % ((bed_res_type,)*6))
#     all_Kmers_ordered = []
#     for (zs, xlinks_z, count_z, occur_xlinks_z, occur_count_z,
#                        xlinks_v, xlinks_v_pval,
#                        count_v, count_v_pval,
#                        occur_xlinks_v, occur_xlinks_v_pval,
#                        occur_count_v, occur_count_v_pval,
#                        xlinks_m, count_m, occur_xlinks_m, occur_count_m,
#                        xlinks_st, count_st, occur_xlinks_st, occur_count_st,
#                        xlinks_cns, count_cns, occur_xlinks_cns, occur_count_cns,
#                        kmer) in o_stats:
# #        xlinks_cns = count_cns = ['intentionally blank']
#         all_Kmers_ordered.append(kmer)
#         fout.write("%s\t%0.3f\t%0.3f\t%0.1f\t%0.3g\t%0.1f\t%0.3g\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%s\t%s\n" % (kmer, occur_xlinks_z, occur_count_z, occur_xlinks_v, occur_xlinks_v_pval, occur_count_v, occur_count_v_pval, occur_xlinks_m, occur_count_m, occur_xlinks_st, occur_count_st, ",".join([str(x) for x in xlinks_cns]), ",".join([str(x) for x in count_cns])))
#     fout.close()
#     print "stats for one (first) occurrences stored to:", out_fname_pref+"_countOne.tab.gz"
#     sys.stdout.flush()
#
#     # write statistics on processed x-links
#     fout = iCount.nc_open(out_fname_pref+"_stats.tab.gz", "wt")
#     fout.write("PARAMETERS\n")
#     fout.write("input file:\t%s\n" % in_fname)
#     fout.write("mapped_to:\t%s\n" % mapped_to)
#     fout.write("annot_ver:\t%s\n" % annot_ver)
#     fout.write("saving to:\t%s\n" % out_fname_pref)
#     fout.write("random permutations:\t%s\n" % random_perms)
#     fout.write("k:\t%s\n" % k)
#     fout.write("intervals:\t%s\n" % intervals)
#     fout.write("report for region:\t%s\n" % str(report_for_region))
#     fout.write("valid region types:\t%s\n" % valid_region_types)
#     fout.write("selected chromes:\t%s\n" % selected_chromes)
#     fout.write("\n")
#
#     fout.write("STATS for [x-link sites]\n")
#     total = sum(xlinks_total.values())
#     fout.write("total:\t%s\n" % total)
#     total = max(total, 1.0)
#     kept = total - sum(xlinks_filtered_out.values())
#     total_filtered = sum(xlinks_filtered_out.values())
#     fout.write("used:\t%s\t%0.1f%%\n" % (kept, kept*100.0/total))
#     fout.write("filtered out:\t%s\t%0.1f%%\t(hit not on same strand as associated gene, not in region(s) selected by user, or no annotation for genomic position)\n" % (total_filtered, total_filtered*100.0/total))
#     fout.write('\n')
#
#     fout.write("region\ttotal\t%\tfiltered\t%\tused\t%\n")
#     for rtype in sorted(set(xlinks_total.keys())|set(xlinks_filtered_out.keys())):
#         reg_total = xlinks_total.get(rtype, 0)
#         reg_filtered = xlinks_filtered_out.get(rtype, 0)
#         reg_kept = reg_total - reg_filtered
#         fout.write("%s\t%s\t%0.1f%%\t%s\t%0.1f%%\t%s\t%0.1f%%\n" % (rtype, reg_total, reg_total*100.0/total, reg_filtered, reg_filtered*100.0/total, reg_kept, reg_kept*100.0/total))
#     fout.write('\n')
#
#     fout.write("STATS for [%s]\n" % bed_res_type)
#     total = sum(count_total.values())
#     fout.write("total:\t%s\n" % total)
#     total = max(total, 1.0)
#     kept = total - sum(count_filtered_out.values())
#     total_filtered = sum(count_filtered_out.values())
#     fout.write("used:\t%s\t%0.1f%%\n" % (kept, kept*100.0/total))
#     fout.write("filtered out:\t%s\t%0.1f%%\t(hit not on same strand as associated gene, not in region(s) selected by users, or no annotation for genomic position)\n" % (total_filtered, total_filtered*100/total))
#     fout.write('\n')
#
#     fout.write("region\ttotal\t%\tfiltered\t%\tused\t%\n")
#     for rtype in sorted(set(count_total.keys())|set(count_filtered_out.keys())):
#         reg_total = count_total.get(rtype, 0)
#         reg_filtered = count_filtered_out.get(rtype, 0)
#         reg_kept = reg_total - reg_filtered
#         fout.write("%s\t%s\t%0.1f%%\t%s\t%0.1f%%\t%s\t%0.1f%%\n" % (rtype, reg_total, reg_total*100.0/total, reg_filtered, reg_filtered*100.0/total, reg_kept, reg_kept*100.0/total))
#     fout.write('\n')
#     fout.close()
#     print "overview stats stored to:", out_fname_pref+"_stats.tab.gz"
#     sys.stdout.flush()
#
#     # distribution of x-linked kmers by gene and by position
#     fout = iCount.nc_open(out_fname_pref+"_by_gene_and_pos.tab.gz", "wt")
#     olst = ["chrome", "s_strand", "start_position", "end_position", "gene_name", "seg_types", "seg_biotypes"]
#     olst.extend(["strand", "xlink_position", bed_res_type] + all_Kmers_ordered)
#     ostr = "\t".join([str(x) for x in olst])
#     fout.write("%s\n" % ostr)
#     # and summary across all positions in each gene
#     sum_by_gene = {}
#     xlinks_by_gene = {}
#     value_by_gene = {}
#     for gene_name in segmentation.genes_order:
#         s_chrome, s_strand, s_segments = segmentation.segments_by_gene[gene_name]
#         if s_strand == '+':
#             start_position = min([sr.min_pos for sr in s_segments])
#             end_position = max([sr.max_pos for sr in s_segments])
#         else:
#             start_position = max([sr.max_pos for sr in s_segments])
#             end_position = min([sr.min_pos for sr in s_segments])
#         seg_types = sorted(list(set(sr.seg_type for sr in s_segments)))
#         seg_types = ",".join(seg_types) if seg_types else "undefined"
#         seg_biotypes = sorted(list(set(sr.biotype for sr in s_segments)))
#         seg_biotypes = ",".join(seg_biotypes) if seg_biotypes else "undefined"
#
#         sum_by_gene[gene_name] = {}
#         tmpd_bypos = bygene_bypos_bykmer.get(gene_name, {})
#         tmpd_bypos = sorted(tmpd_bypos.items())
#         for (chrome, strand, p_start, value, seg), tmpd_bykmer in tmpd_bypos:
#             xlinks_by_gene[gene_name] = xlinks_by_gene.get(gene_name, 0) + 1
#             value_by_gene[gene_name] = value_by_gene.get(gene_name, 0) + value
#             assert(chrome == s_chrome)
#             for kmer, cn in tmpd_bykmer.iteritems():
#                 sum_by_gene[gene_name][kmer] = sum_by_gene[gene_name].get(kmer, 0) + cn
#             olst = [chrome, s_strand, str(start_position), str(end_position), gene_name, seg_types, seg_biotypes]
#             olst.extend([strand, p_start, value] + [str(tmpd_bykmer.get(kmer, '')) for kmer in all_Kmers_ordered])
#             ostr = "\t".join([str(x) for x in olst])
#             fout.write("%s\n" % ostr)
#     fout.close()
#
#     # summary by gene
#     fout = iCount.nc_open(out_fname_pref+"_by_gene_summary.tab.gz", "wt")
#     olst = ["chrome", "s_strand", "start_position", "end_position", "gene_length", "gene_name", "seg_types", "seg_biotypes", "xlinks", "xlinks density", "%s sum" % bed_res_type, "%s sum density" % bed_res_type]
#     for kmer in all_Kmers_ordered:
#         olst.append("%s freq." % kmer)
#         olst.append("%s density" % kmer)
#     ostr = "\t".join([str(x) for x in olst])
#     fout.write("%s\n" % ostr)
#     for gene_name in segmentation.genes_order:
#         s_chrome, s_strand, s_segments = segmentation.segments_by_gene[gene_name]
#         if s_strand == '+':
#             start_position = min([sr.min_pos for sr in s_segments])
#             end_position = max([sr.max_pos for sr in s_segments])
#         else:
#             start_position = max([sr.max_pos for sr in s_segments])
#             end_position = min([sr.min_pos for sr in s_segments])
#         gene_len = float(abs(end_position - start_position) + 1)
#         seg_types = sorted(list(set(sr.seg_type for sr in s_segments)))
#         seg_types = ",".join(seg_types) if seg_types else "undefined"
#         seg_biotypes = sorted(list(set(sr.biotype for sr in s_segments)))
#         seg_biotypes = ",".join(seg_biotypes) if seg_biotypes else "undefined"
#
#         tmpd_bykmer = sum_by_gene.get(gene_name, {})
#         olst = [s_chrome, s_strand, str(start_position), str(end_position), str(gene_len), gene_name, seg_types, seg_biotypes]
#         olst.extend([str(xlinks_by_gene.get(gene_name, 0)), "%0.5f" % (xlinks_by_gene.get(gene_name, 0)/gene_len)])
#         olst.extend([str(value_by_gene.get(gene_name, 0)), "%0.5f" % (value_by_gene.get(gene_name, 0)/gene_len)])
#         for kmer in all_Kmers_ordered:
#             if kmer in tmpd_bykmer:
#                 olst.append(tmpd_bykmer[kmer])
#                 olst.append("%0.5f" % (tmpd_bykmer[kmer]/gene_len))
#             else:
#                 olst.append('')
#                 olst.append('')
#         ostr = "\t".join([str(x) for x in olst])
#         fout.write("%s\n" % ostr)
#     fout.close()
#
#     if analysis_id is not None:
#         files = [
#             ('tab-stats',
#              iCount.strip_storage_prefix(out_fname_pref+"_stats.tab.gz"),
#              db.filesize_nicefmt(out_fname_pref+"_stats.tab.gz")
#             ),
#             ('tab-allOccurrences',
#              iCount.strip_storage_prefix(out_fname_pref+"_countAll.tab.gz"),
#              db.filesize_nicefmt(out_fname_pref+"_countAll.tab.gz")
#             ),
#             ('tab-oneOccurrence',
#              iCount.strip_storage_prefix(out_fname_pref+"_countOne.tab.gz"),
#              db.filesize_nicefmt(out_fname_pref+"_countOne.tab.gz")
#             ),
#             ('tab-byGeneAndPos',
#              iCount.strip_storage_prefix(out_fname_pref+"_by_gene_and_pos.tab.gz"),
#              db.filesize_nicefmt(out_fname_pref+"_by_gene_and_pos.tab.gz")
#             ),
#             ('tab-byGeneSummary',
#              iCount.strip_storage_prefix(out_fname_pref+"_by_gene_summary.tab.gz"),
#              db.filesize_nicefmt(out_fname_pref+"_by_gene_summary.tab.gz")
#             ),
#         ]
#         db.analysis_set_output(analysis_id, files)
#         db.analysis_set_status(analysis_id, 'done')
#         db.analysis_set_last_update(analysis_id)
#     print "DONE"
#     print time.ctime()
#     print "*********************************"
#     print
#     sys.stdout.flush()
#
