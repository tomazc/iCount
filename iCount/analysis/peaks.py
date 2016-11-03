"""
Peak finding
------------

Find positions with high density of cross-linked sites.

There are two typical variants of this analysis, depending on the parameters:

    * Gene-wise analysis, where:
        * features = gene
        * group_by = gene_id
    * Transcript-wise analysis where:
        * features = CDS, intron, UTR3, UTR5, ncRNA, intergenic
        * group_by = transcript_id


Let's look at the Gene-wise analysis in more detail first. Imagine the following
situation::

    |-----------gene1----------|
            |-----------------------------gene2------------------------------|
           ab c                d                 e

    a = 60
    b = 100
    c = 70
    d = 40
    e = 100
    gene1: gene_id = 001
    gene2: gene_id = 002

There are two genes (partially intersecting) and five positions with cross-links
(noted with a, b, c, d and e). Crosslink position "a" has 60 cross-link events,
"b" has 100 cross-link events and so on. Also, gene1 has gene_id 001, etc.

The algorithm first finds all intersections between annotation and cross-links.
In this case cross-link position "a" intersects only with gene1, while position
"b" intersects also with gene2... Annotation can include various other types of
segments (transcripts, intergenic, ncRNA, etc.), but only segments of type
``gene`` are considered for intersection. This behaviour is controlled by
parameter ``features``.

Next step is to make groups of cross-links. They are grouped by ``group_by``
parameter (in this case, it equals to ``gene_id``). There will be 2 groups.
First group name will be 001 and will contain a, b, c and d. Second group name
will be 002 and will contain b, c, d and e.

The question now is: has any of positions in each group significantly increased
number of cross-link events? And how can one quantify this significance?

This is done by parmutation analysis. It draws a number of random situations
with same group size and number of cross-link scores. Number of such draws is
determined by ``perm`` parameter. This way, a random distribution is calculated.
When comparing the observed distribution with the random one, FDR values are
assigned to each position. A cutoff FRD value is chosen and only positions with
FDR < FDR cutoff are considered as significant.

One must also know that when considering only scores on single positions
significant *clusters* of cross-links can be missed. In the upper example, it is
obviuous, that something more significantly is happening on position b than on
position e, despite having the same score. To account for this, algorithm
considers not only the score of single cross-link, but also scores of
cross-links some nucleotides before and after. This behaviour in controlled by
half-window (hw) parameter. In the upper example, score of position b eqals to
160 if hw = 1 and 2530 if hw=2. Score of position e remains 100.

Let's also look at the transcript-wise analysis. In this case, scenario also
includes transcripts and sub-transcript elements::


    |-----------gene1----------|
    |--------transcript1-------|
    |-ncRNA-||-intron-||-ncRNA-|
            |-----------------------------gene2------------------------------|
            |---------------transcript2--------------|
            |-CDS-||-intron-||-CDS-||-intron-||-UTR3-|
                                    |---------------transcript3--------------|
                                    |-UTR5-||-intron-||-CDS-||-intron-||-CDS-|

           ab c                d                 e

    a = 60
    b = 100
    c = 70
    d = 40
    e = 100
    gene1: gene_id = 001
    gene2: gene_id = 002
    transcript1: transcript_id = 0001
    transcript2: transcript_id = 0002
    transcript3: transcript_id = 0003

Value of parameter features is: CDS, intron, UTR3, UTR5, ncRNA, intergenic.
Value of parameter group_by is transcript_id. Since we have multiple values in
feature parameter, another parameter becomes important: merge_features. If set
to false (default) algorithm will make the following groups:

    * group name: ncRNA-0001, members: a, b, d
    * group name: intron-0001, members: c
    * group name: CDS-0002, members: b, c, d
    * group name: UTR3-0002, members: e
    * group name: intron-0003, members: e

However, if merge_features equals to true, groups are:

    * group name: 0001, members: a, b, c, d
    * group name: 0002, members: b, c, d, e
    * group name: 0003, members: e

Then, for each group, procedure is exactly the same as in Gene-wise case.

When analysis is done, significant positions are reported in file, given by
peaks parameter. If scores parameter is also given, all positions are
reported in it, no matter the FDR value.

"""
import math
import bisect
import logging

import numpy
import pybedtools

import iCount

from collections import Counter

from iCount.files.bed import _f2s

LOGGER = logging.getLogger(__name__)


def _sum_within_window(pos_val, hw=3):
    """
    Sum counts in windows of half-window size ``hw`` in ``pos_val``.

    Example::

        pos_val = [(41, 1), (42, 1), (43, 1), (44, 1)]
        ret_list = _sum_within_window(pos_val, 1)
        pos_val = [(41, 2), (42, 3), (43, 3), (44, 2)]

        pos_val = [(41, 1), (42, 1), (43, 1), (46, 1)]
        ret_list = _sum_within_window(pos_val, 3)
        pos_val = [(41, 3), (42, 3), (43, 4), (44, 2)]

    The returned list preserves the order of positions given on input.
    """
    if not pos_val:
        return []
    pos_val_ind = sorted((p, v, i) for i, (p, v) in enumerate(pos_val))
    poss, vals, inds = zip(*pos_val_ind)

    ret_list = [None] * len(pos_val)
    max_i = len(poss)
    for i, p in enumerate(poss):
        i_start = bisect.bisect_left(poss, p - hw, lo=max(0, i - hw), hi=i)
        i_stop = bisect.bisect_left(poss, p + 1 + hw, lo=i_start, hi=min(i + hw + 1, max_i))
        ret_list[inds[i]] = (p, sum(vals[i_start:i_stop]))
    return ret_list


def _sum_within_window_nopos(pos_val, hw=3):
    """
    Same as _sum_within_window but without positions.
    """
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
    """
    Given a list of SWW scores in region ``vals``, return list ``freqs_cum`` where
    probability that randomly picked value in `vals` is equal or greater than i
    equals freqs_cum[i].

    Max_val is the largest possible value that can be expected in `vals`.
    """
    # Make histogram from vals, with max_val + 2 bins. O bin for zero scores and
    # one for one more than max_val. The last one bin shoud alwayxs be zero, right?

    # + 1 for the probability of observing 0 or more
    # + 1 beacouse range() excludes last element
    freqs, _ = numpy.histogram(vals, bins=range(max_val+1+1), density=True)

    # Now we want to know not how many events with exactly x cross links is
    # possible, but with x cross-links OR MORE. We sum from behind:
    freqs_cum = numpy.cumsum(freqs[::-1])

    # Still, we return the reversed freqs_cum.
    return freqs_cum[::-1]


# def get_rnd_distrib(size, total_hits, hw, perms=100):
#     """The simplest (and fastest) permutation test.
#
#     Not used by default.
#     """
#     rnd_cns = numpy.zeros(total_hits+1)
#     for i in range(perms):
#         rnd_hits = Counter(numpy.random.randint(size, size=total_hits))
#         rnd_hits_extended = _sum_within_window_nopos(rnd_hits.items(), hw)
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


def get_avg_rnd_distrib(size, total_hits, hw, perms=10000):
    """
    Return background distribution of peak heights for given region size
    and number of hits.

    We follow the modified FDR for peak height, proposed by [1]

    [1] Yeo, G.W. et al. An RNA code for the FOX2 splicing regulator revealed
    by mapping RNA-protein interactions in stem cells. Nat. Struct. Mol.
    Biol. 16, 130–137 (2009).
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2735254/

    Results are cached, so they can be reused.

    Parameters
    ----------
    size : int
        Size of region.
    total_hits : int
        Number of cross-link events in region.
    hw : int
        Half-window size. The actal window size is: 2 * hw + 1.
    perms : int
        Number of permutations to make.

    Returns
    -------
    numpy.ndarray
        Probability to find CWW score i or more on chosen position is equal to
        i-th element or returned array.
    """
    global ps_cache
    cache_key = (size, total_hits, hw, perms)
    if cache_key not in ps_cache:

        rnd_ps = numpy.zeros((perms, total_hits + 1))
        for i in range(perms):

            # Draw random distribution of cross-link events in a group with
            # group size = `size` and number of cross-link events = `total_hits`
            rnd_hits = Counter(numpy.random.randint(size, size=total_hits))

            # This is then list. i-th element in list is probability, that there
            # is equal or more than i crossslinks on some position???

            scores_cww = _sum_within_window_nopos(rnd_hits.items(), hw=hw)
            rnd_ps[i, :] = cumulative_prob(scores_cww, total_hits)

        rnd_dist = numpy.mean(rnd_ps, axis=0) + numpy.std(rnd_ps, axis=0)
        # Adding std, can make probability higher than 1, which is nonsense. Fix:
        rnd_dist_fixed = [min(1.0, prob) for prob in rnd_dist]
        ps_cache[cache_key] = rnd_dist_fixed

    return ps_cache[cache_key]


def _process_group(pos_scores, group_size, hw, perms):
    """
    Assign FDR value to each position in group.

    One is given a region of size N, with given distribution of K cross-link
    events (scores). How likely it is, that the observed distribution of scores
    is random? Additionally, scores need to be averaged in windows of size `hw`.

    Lets do an example, with region size 5 and 3 cross-link events::

         0 0 2 0 1
        |_|_|_|_|_|

        pos_scores = [(0, 0), (1, 0), (2, 2), (3, 0), (4, 1)]

    First, do the averaging. This is done by function ``_sum_within_window``. If
    half-window size is hw=1, the result is::

         0 2 2 3 1
        |_|_|_|_|_|

        pos_scores_sww = [(0, 0), (1, 2), (2, 2), (3, 3), (4, 1)]

    For *this* particular situation, what is the probability, to find a position
    in region with 3 or more SWW scores? Well, only one position has SWW
    score 3 or more. Since we have 5 positions, the probability of observing SWW
    score of 3 or more is 1/5=0.2. What about 2 or more? We have 2 positions
    with SWW score 2 and one position with score 3. So, in total we have 3
    positions. Probability of observing position with 2 or more SWW scores is
    3 / 5 = 0.6. This way, one can make the calculation for all cases::

        Prob (>= 0) = 5 / 5 = 1.0
        Prob (>= 1) = 4 / 5 = 0.8
        Prob (>= 2) = 3 / 5 = 0.6
        Prob (>= 3) = 1 / 5 = 0.2

    These probabilities can be put in a list: [1.0, 0.8, 0.6, 0.2]. Storing
    results in such list has a nice property::

        Prob (>= i) = list[i]

    ``cumulative_prob`` function produces exactly such list from SWW scores.
    This is "cumulative_prob" for given example.

    To produce a reference to which this can be compared, we use function
    ``get_avg_rnd_distrib``.

    Then, random and observed "cumulative_prob" are compared and FDR scores for
    each cross-link can be derived. More can be read in artice [1] or in code
    comments.

    [1] Yeo, G.W. et al. An RNA code for the FOX2 splicing regulator revealed by
    mapping RNA-protein interactions in stem cells. Nat. Struct. Mol. Biol. 16,
    130–137 (2009).

    Parameters
    ----------
    pos_scores : list
        Lits with (position, scores) elements.
    group_size : list
        Size of region

    Returns
    -------
    list
        List of tuples, containing (position, score, sww_score and fdr_value)
        for all cross-link positions in group.

    """
    # Count the number of all cross-link events in a group:
    sum_scores = math.ceil(sum([score for _, score in pos_scores]))

    # Calculate the observed cumulative_prob:
    pos_scores_sww = _sum_within_window(pos_scores, hw=hw)
    positions, scores_sww = zip(*pos_scores_sww)
    observed = cumulative_prob(scores_sww, sum_scores)

    # Calculate random cumulative_prob for given group_size and sum_scores:
    random_ = get_avg_rnd_distrib(group_size, sum_scores, hw, perms=perms)

    # This step follows the article [1] to produce FDR values. First, produce
    # mapping from sww_scores to FDR value:
    sww2fdr = [min(rnd / true, 1.0) for rnd, true in zip(random_, observed)]
    # Compute FDR scores por each position based on it's sww_score:
    fdr_scores = [sww2fdr[round(sww_score)] for _, sww_score in pos_scores_sww]

    positions, scores = zip(*pos_scores)
    return zip(positions, scores, scores_sww, fdr_scores)


def run(annotation, sites, peaks, scores=None, features=None, group_by='gene_id',
        merge_features=False, hw=3, fdr=0.05, perms=100, rnd_seed=42, report_progress=False):
    """
    Find positions with high density of cross-linked sites.

    Algorithm does:

        * read the annotation file
        * read the BED file with cross-links



    When detrmining feature.name, value of the first existing attribute in the
    following tuple is taken::

        ("ID", "gene_name", "transcript_id", "gene_id", "Parent")

    Source in pybedtools:
    https://github.com/daler/pybedtools/blob/master/pybedtools/scripts/annotate.py#L34

    Parameters
    ----------
    annotation : str
        Annotation file in GTF format
    sites : str
        File with cross-links in BED6 format.
    peaks : str
        File name for "peaks" output. File reports positions with significant
        number of cross-link events. It should have .bed or .bed.gz extension.
    scores : str
        File name for "scores" output. File reports all cross-link events,
        independent from their FDR score It should have .tsv, .csv, .txt or .gz
        extension.
    features : list_str
        Features from annotation to consider. If None, 'gene' is used.
    group_by : str
        Attribute by which cross-link positions are grouped.
    merge_features : bool
        Treat all features as one when grouping. Has no effect when only one
        feature is given in features parameter.
    hw : int
        Half-window size.
    fdr : float
        FDR threshold.
    perms : int
        Number of permutations when calculating random distribution.
    rnd_seed : int
        Seed for random generator.
    report_progress : bool
        Report analysis progress.

    Returns
    -------
    iCount.metrics
        Analysis metadata.

    """
    iCount.log_inputs(LOGGER, level=logging.INFO)

    if fetaures is None:
        fetaures = ['gene']
    assert peaks.endswith(('.bed', '.bed.gz'))
    if scores:
        assert scores.endswith(('.tsv', '.tsv.gz', '.csv', '.csv.gz', 'txt', 'txt.gz'))
    numpy.random.seed(rnd_seed)

    LOGGER.info('Loading annotation file...')
    annotation = pybedtools.BedTool(annotation)
    LOGGER.info('Loading cross-links file...')
    sites = pybedtools.BedTool(sites).sort().saveas()

    # intersect cross-linked sites with regions
    LOGGER.info('Calculating intersection between annotation and cross-link file...')
    overlaps = annotation.intersect(sites, sorted=True, s=True, wo=True).saveas()

    groups = {}
    group_sizes = {}
    metrics = iCount.Metrics()
    metrics.skipped_features = 0  # conuter for skipped fetaures
    multi_mode = len(features) > 1 and not merge_features
    LOGGER.info('Processing intersections...')
    for feature in overlaps:
        if feature[2] not in features:  # check if of correct feature type
            metrics.skipped_features += 1
            continue

        chrom = feature.chrom
        start = feature.start
        end = feature.stop
        name = feature.name
        strand = feature.strand
        site_pos = int(feature.fields[10])
        site_score = float(feature.fields[13])
        # Determine group_id depending on multi_mode...
        group_id = feature.attrs[group_by]
        if multi_mode:
            group_id = feature[2] + '_' + group_id

        groups.setdefault((chrom, strand, group_id, name), []).append((site_pos, site_score))
        group_sizes.setdefault((chrom, strand, group_id, name), set()).add((start, end))

    # Validate that segments in same group do not overlap: start of next feature
    # is greater than stop of the current one:
    for sizes in group_sizes.values():
        sizes = sorted(sizes)
        for first, second in zip(sizes, sizes[1:]):
            assert first[1] < second[0]

    # calculate total length of each group by summing element sizes:
    group_sizes = dict([(name, sum([end - start for start, end in elements])) for
                       name, elements in group_sizes.items()])

    # calculate and assign FDRs to each cross-linked site. FDR values are
    # calcualated together for each group.
    results = {}
    metrics.all_groups = len(groups)
    progress, j = 0, 0
    for (chrom, strand, group_id, name), hits in sorted(groups.items()):
        j += 1
        if report_progress:
            new_progress = j / metrics.all_groups
            progress = iCount._log_progress(new_progress, progress, LOGGER)

        group_size = group_sizes[(chrom, strand, group_id, name)]

        # Crucial step: each position in a group is given a fdr_score, based on
        # hits in group, group_size, half-window size and number of
        # permutations. Than, FDR scores (+ some other info) are written to
        # `results` container:
        processed = _process_group(hits, group_size, hw, perms)
        for (pos, val, val_extended, fdr_score) in processed:
            results.setdefault((chrom, pos, strand), []).\
                append((fdr_score, name, group_id, val, val_extended))

    LOGGER.info('Peaks caclulation finished. Writing results to files...')

    # Make peaks: a BED6 file, with only the most significant cross-links:
    metrics.significant_positions = 0
    with iCount.files.gz_open(peaks, 'wt') as peaks:
        for (chrom, pos, strand), annot_list in sorted(results.items()):
            annot_list = sorted(annot_list)

            # report minimum fdr_score for each position in BED6
            min_fdr_score = annot_list[0][0]
            if min_fdr_score < fdr:
                metrics.significant_positions += 1
                # position has significant records - report the most significant ones:
                min_fdr_records = [rec for rec in annot_list if rec[0] == min_fdr_score]

                _, names, group_ids, group_scores, _ = zip(*min_fdr_records)
                name = ','.join(names) + ' - ' + ','.join(group_ids)
                line = [chrom, pos, pos + 1, name, group_scores[0], strand]
                peaks.write('\t'.join([_f2s(i, dec=4)for i in line]) + '\n')
    LOGGER.info('BED6 file with significant peaks saved to %s', peaks)

    # Make scores: a tab-separated file, with ALL cross-links, (no significance threshold)
    header = ['chrom', 'position', 'strand', 'name', 'group_id', 'score', 'score_extended', 'FDR']
    if scores:
        with iCount.files.gz_open(scores, 'wt') as scores:
            scores.write('\t'.join(header) + '\n')
            for (chrom, pos, strand), annot_list in sorted(results.items()):
                for (fdr_score, name, group_id, score, val_extended) in sorted(annot_list):
                    line = [chrom, pos, strand, name, group_id, score, val_extended, fdr_score]
                    scores.write('\t'.join([_f2s(i, dec=6) for i in line]) + '\n')
        LOGGER.info('Scores for each cross-linked position saved to %s', scores)

    LOGGER.info('Done.')
    return metrics
