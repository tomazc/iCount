""".. Line to protect from pydocstyle D205, D400.

RNA maps
--------

Distribution of cross-links relative to genomic landmarks.

What is an RNA-map?
^^^^^^^^^^^^^^^^^^^

Imagine the following situation (on positive strand)::

    |---intron1---||---exon2---||---intron2---||---exon3---|
                              x|-----R1-----|

Situation is simple: a read perfectly maps to reference genome. First cross-link
(nucletide before the start of read R1) is located -2 nucleotides before
landmark of type 'exon-intron'. By examining all cross-links(reads) that are
located in vicinity of landmarks one can obtain the distribution of
cross-link-to-landmark distances. Such distibution is called an RNA map.

Of course, situation gets complicated. Real-world-situation more likely looks
like this::

    |--------transcript1-------|
    |-ncRNA-||-intron-||-ncRNA-|
            |------------------transcript2-----------------|
            |--CDS--||-intron-||--CDS--||-intron-||--UTR3--|
                                    |---------------transcript3--------------|
                                    |-UTR5-||-intron-||-CDS-||-intron-||-CDS-|

             x|-----R1-----|
             x|R2.1->          <--R2.2-|
                                x|-R3-|


All sort of situations can arise:

    * read can intersect with multiple genes/transcripts/segments
      simultaneously, which means that one read is relevant for different types
      of RNA maps.
    * if whole read is mapped to one segment (as is read R3 on transcript2),
      then it is impossible to tell the type of RNA map it belongs to: is it
      CDS-intron or intron-CDS?
    * read can be mapped in multiple parts - for example when introns are
      removed before translation, first half of read will map to ``exon1`` and
      second half to ``exon2``.
    * same cross-link event can be represented by multiple reads with same
      random barcode.
    * ...


The algorithm takes care for all of this. Finally, it outputs a container
``data`` with the following structure::

    data = {
        rna_map_type: {relative_position: [all, explicit], ...}
        exon-intron: {-42: [1042, 114], -41: [784, 69] ...}
        exon-intron-neg: {-42: [219, 65], ...}
        ...
    }

The structure is a dict, with keys corresponding to different RNA map types.
Values are also dicts. Keys in these dicts are positions of cross links relative
to the landmark of the RNA map. Values for each position are lists: the first
and second element count number of *all* and *explicit* cross-link events
respectively. Cross link is explicit, if read identifying cross-link is mapped
to both parts of RNA-map. In upper example, R1 is explicit, since it maps to
intron *and* exon. Read R3 on transcript2 is implicit, since one has to decide
wheather it's cross-link belongs to ``CDS-intron`` or ``intron-CDS`` RNA map.
The term *all* means explicit + implicit in this context.

If read is implicit (start and stop within same segment), one can choose two
varinats of the algorithm. One can choose that whole score of read is given to
RNA-map of the closest landmark. Other option is to split score on half to both
neighbouring segments. This behaviour is controlled by parameter
implicit_handling.

There are also two cases when read can map in a way that is not predicted with annotation:

    * For split reads, second start can fall on nucleotide that is not start of
      a new segment. This behaviour is excactly the same as in xlsites analysis.
      Such reads are reported in file defined with ``strange`` parameter.
    * Reads that map on two different transcripts or even two different genes
      are repoted in file with name defined in cross_transcript.

-------------------------------------------------------------------------------

Things that still need to be resolved:

    * negative strand inspection:
        * for each ss_group, take annotation from reverse strand and check what
          is the situation on other side.
        * What si again the bio-contxt here - explain in docs...?
    * Test: Each cross-link should have the same cDNA count as in xlsites
      results. The final check should be: sum all scores - the result should be
      the same as for xlsites. Actually, the metrics.cross_transcript should
      also be considered:
      assert: sum(xlsites scores) == sum(RNAmaps "all" + metrics.cross_transc)
    * if read crosses two ore more regions from start to end, only end-type is
      considered for RNAmap type. Should we fix this and which region type to
      take? Since reads are tpyically much shorted than regions, this is
      probably a rare case.

-------------------------------------------------------------------------------

Wishes and ideas(Jernej, Tomaz):

    * The whole idea is to get a feeling what is happening around landmarks

    * separate pre-mRNA and mRNA [this is partially done in _add_entry]
        * Life cyle of RNA: part of DNA in transcribed to RNA in nucleus. Introns
          are cut out already between this process! Really???
        * pre-mRNA contains introns, mature mRNA has no introns and goes outside
          from the nucleus.
        * so: RNAmap exon-intron will for sure contain only pre-mRNA, but RNAmap
          exon-exon can contain pre-mRNA and mRNA... It would be really nice to
          report on percentage of theese two categories.
        * This can be done by introducing three categories:  pre-mRNa, mRNA and
          undecided. This can be determined from results - they contain explicit /
          implicit info and RNAmap type:

    * Intoduce three categories (colors in visualisation?) in every RNAmap type.
      Let's take exon-intron as example:

        * whole read in left part - color #1 (negative coord, implicit)
        * whole read in right part - color #2 (positive coord, implicit)
        * read crossing junction - color #3 (explicit)

    * How to handle situation when one region in individual's sequence clearly
      differs from reference sequence but it is just some variation?

        * Change the reference sequence? This can be complex... Make a helper
          tool for this?
        * Provide additional data in function - which exceptions /abnormalities
          to ignore?

    * Gaussian smooting of RNAmaps?

    * split read - it can differ also on "first-stop", not juts second-start
        * we could have some sort of quality check when observing the variance
          of first-stop-s. THe major question is: "Do all reads for certain
          xlink event indicate the same behaviour?"

    * Normalization?

"""
import logging

import iCount
from iCount.files import _f2s

LOGGER = logging.getLogger(__name__)

EXON_TYPES = ['CDS', 'ncRNA', 'UTR3', 'UTR5']


def _add_entry(start_type, stop_type, dist, score, strand, data, metrics, explict=False):
    """Add RNA-map entry in ``data``."""
    if strand == '+':
        rna_map_type = start_type + '-' + stop_type
    else:
        rna_map_type = stop_type + '-' + start_type
        dist = -dist

    data.setdefault(rna_map_type, {}).setdefault(dist, [0, 0])[0] += score
    if explict:
        data.setdefault(rna_map_type, {}).setdefault(dist, [0, 0])[1] += score

    # Increment also mRNA / pre-mRNa counters:
    if 'intron' in rna_map_type or 'intergenic' in rna_map_type:
        metrics.origin_premrna += score
    elif start_type in EXON_TYPES and stop_type in EXON_TYPES:
        metrics.origin_mrna += score
    else:
        metrics.origin_ambiguous += score


def _process_read_group(xlink, chrom, strand, read, data, annotation, metrics,
                        implicit_handling='closest'):
    """
    Process each read group.

    Each group is represented by read. It has a fixed xlink, start, second-start
    and end position. This function finds all intersections with annotation and
    determines for each of them:

        * type of RNA map (+ explicit/implicit)
        * relative position to landmark
        * score that belongs to each intersection


    Note#1: ``annotation`` contains just the genes that span over cross-link
    position (+ one before and one after). It has the following shape::

        annotation = {
            gene_id#1: {
                'gene_segment': gene_segment,
                transcript_id#1: [transcript_segment, exon1, intron1, exon2, ...],
                transcript_id#2: [transcript_segment, exon1, intron1, exon2, ...],
                ...
            },
            gene_id#2: {},
            ...
        },

    Note#2: ``read`` is a tuple with the following content::

        read_data = (middle_pos, end_pos, read_len, num_mapped, second_start)

    """
    stop = read[1]  # stop is in second column
    start = xlink + (1 if strand == '+' else - 1)

    containing_start = {}
    containing_stop = {}
    for _, gene_content in annotation:
        for transcript_id, transcript_content in gene_content.items():
            if transcript_id == 'gene_segment':
                continue

            transcript_segment = transcript_content[0]
            if transcript_segment.start <= start <= transcript_segment.stop:
                containing_start[transcript_id] = transcript_content
            if transcript_segment.start <= stop <= transcript_segment.stop:
                containing_stop[transcript_id] = transcript_content

    containing_both = set(containing_start.keys()) & set(containing_stop.keys())
    # There are transcripts that contin start and stop of the read:
    if containing_both:
        relevant_transcripts = dict(
            [trs for _, g_cnt in annotation for trs in g_cnt.items() if trs[0] in containing_both])
    # There are NO transcripts that contin start and stop of the read. Remaining
    # options that do not violate annotation are:
    # (a) 'intergenic-transcript' RNA map type:
    elif len(containing_start) == 1 and list(containing_start.values())[0][0][2] == 'intergenic':
        tr_score = 1 / len(containing_stop)
        for transcript_id, transcript_content in containing_stop.items():
            stop_segment_index = next(
                (i for i, segment in enumerate(transcript_content) if
                 segment.start <= stop <= segment.stop and segment[2] != 'transcript'))
            start_type = 'intergenic'
            stop_type = transcript_content[stop_segment_index][2]
            rel_dist = xlink - transcript_content[stop_segment_index].start
            _add_entry(
                start_type, stop_type, rel_dist, tr_score, strand, data, metrics, explict=True)

        return  # skip the rest of the algorithm

    # (b) 'transcript-intergenic' RNA map type:
    elif len(containing_stop) == 1 and list(containing_stop.values())[0][0][2] == 'intergenic':
        tr_score = 1 / len(containing_start)
        for transcript_id, transcript_content in containing_start.items():
            start_segment_index = next(
                (i for i, segment in enumerate(transcript_content) if
                 segment.start <= start <= segment.stop and segment[2] != 'transcript'))
            start_type = transcript_content[start_segment_index][2]
            stop_type = 'intergenic'
            rel_dist = xlink - transcript_content[start_segment_index].stop
            _add_entry(
                start_type, stop_type, rel_dist, tr_score, strand, data, metrics, explict=True)

        return  # skip the rest of the algorithm

    # This read has mapped in way that is not predicted by annotation: it should be reported:
    else:
        # TODO: Ideally, this would produce a BAM file with such reads... but
        # read instance from all BAM related info is far back in
        # _processs_bam_file function... For now:
        metrics.cross_transcript += 1
        data.setdefault('cross_transcript', {}). \
            setdefault((chrom, strand, xlink), []).append(read)
        return

    # ###################################################

    # Note that only "containing_both" scenario reaches this point.

    tr_score = 1 / len(relevant_transcripts)
    relevant_transcripts = sorted(relevant_transcripts.items(), key=lambda x: x[1][0].start)
    for transcript_id, transcript_content in relevant_transcripts:
        transcript_content = [s for s in transcript_content if s[2] != 'transcript']
        start_segment_index = next((i for i, segment in enumerate(transcript_content)
                                    if segment.start <= start <= segment.stop))
        stop_segment_index = next((i for i, segment in enumerate(transcript_content)
                                   if segment.start <= stop <= segment.stop))

        # Explicit case: this is easy
        if start_segment_index != stop_segment_index:
            start_type = transcript_content[start_segment_index][2]
            stop_type = transcript_content[stop_segment_index][2]
            rel_dist = xlink - transcript_content[stop_segment_index].start
            _add_entry(
                start_type, stop_type, rel_dist, tr_score, strand, data, metrics, explict=True)

        # Implicit case: this can be tricky...
        else:
            # Container for all options of RNA map type:
            options = []

            segment = transcript_content[start_segment_index]
            rel_dist_down = xlink - segment.start
            rel_dist_up = xlink - segment.stop
            # Note: segment_n.stop == segment_n+1.start

            # ###################################################

            # Handle downstream
            if start_segment_index != 0:
                start_type = transcript_content[start_segment_index - 1][2]
                options.append([start_type, segment[2], rel_dist_down])
            else:
                # Gene beefore start (downstream) is the first entry in annotation:
                gene_down = annotation[0][1]['gene_segment']
                # this segment OR the downstream gene has to be intergenic for
                # this to be OK with annotation:
                if 'intergenic' in [gene_down[2], segment[2]] and gene_down.stop == segment.start:
                    trs_down = [tr_cnt for tr_id, tr_cnt in annotation[0][1].items() if
                                tr_id != 'gene_segment' and tr_cnt[0].stop == segment.start]
                    for tr_cnt in trs_down:
                        seg_down = next((seg for seg in tr_cnt if seg.stop == segment.start and
                                         seg[2] != 'transcript'))
                        options.append([seg_down[2], segment[2], rel_dist_down])
                else:
                    # Transcript-transscript scenario - not allowed. (no appends to options)
                    pass

            # ###################################################

            # Handle upstream: (similar to downstream)
            if stop_segment_index != len(transcript_content) - 1:
                stop_type = transcript_content[start_segment_index + 1][2]
                options.append([segment[2], stop_type, rel_dist_up])
            else:
                gene_up = annotation[-1][1]['gene_segment']
                if 'intergenic' in [gene_up[2], segment[2]] and gene_up.start == segment.stop:
                    trs_up = [tr_cnt for tr_id, tr_cnt in annotation[-1][1].items() if
                              tr_id != 'gene_segment' and tr_cnt[0].start == segment.stop]
                    for tr_cnt in trs_up:
                        seg_up = next((seg for seg in tr_cnt if seg.start == segment.stop and
                                       seg[2] != 'transcript'))
                        options.append([segment[2], seg_up[2], rel_dist_up])
                else:
                    pass

            # ###################################################

            # Handle also exons (if segment type is exon and introns are removed
            # there is possibility that rna map is also of exon-exon type):
            if segment[2] in EXON_TYPES:
                exon_number = int(segment.attrs['exon_number'])
                exon_before = next((s for s in transcript_content if 'exon_number' in s.attrs and
                                    int(s.attrs['exon_number']) == exon_number - 1), None)
                exon_after = next((s for s in transcript_content if 'exon_number' in s.attrs and
                                   int(s.attrs['exon_number']) == exon_number + 1), None)
                if exon_before:
                    options.append([exon_before[2], segment[2], rel_dist_down])
                if exon_after:
                    options.append([segment[2], exon_after[2], rel_dist_up])

            # ###################################################

            if implicit_handling == 'closest':
                # Compute minimal distance from border:
                min_dist = options[0][2]
                for _, _, dist in options:
                    if abs(dist) < abs(min_dist):
                        min_dist = dist

                # Take all options that are on minimal absolute distance:
                options = [i for i in options if i[2] == min_dist]

            # Write entries in options to ``data``:
            final_score = tr_score / len(options)
            for start_type, stop_type, rel_dist in options:
                _add_entry(start_type, stop_type, rel_dist, final_score, strand, data, metrics)


def run(bam, segmentation, out_file, strange, cross_transcript, implicit_handling='closest',
        mismatches=2, mapq_th=0, holesize_th=4):
    """
    Compute distribution of cross-links relative to genomic landmarks.

    Parameters
    ----------
    bam : str
        BAM file with alligned reads.
    segmentation : str
        GTF file with annotation. Should be a file produced by function
        `get_regions`.
    out_file : str
        Output file with analysis results.
    strange : str
        File with strange propertieas obtained when processing bam file.
    cross_transcript : str
        File with reads spanning over multiple transcripts or multiple genes.
    implicit_handling : str
        Can be 'closest' or 'split'. In case of implicit read - split score to
        both neighbours or give it just to the closest neighbour.
    mismatches : int
        Reads on same position with random barcode differing less than
        ``mismatches`` are grouped together.
    mapq_th : int
        Ignore hits with MAPQ < mapq_th.
    holesize_th : int
        Raeads with size of holes less than holesize_th are treted as if they
        would have no holes.

    Returns
    -------
    str
        File with number of (al, explicit) scores per each position in each
        RNA-map type.

    """
    iCount.logger.log_inputs(LOGGER)

    if implicit_handling not in ('closest', 'split'):
        raise ValueError(
            'Parameter implicit_handling should be one of "closest" or "split"')

    metrics = iCount.Metrics()
    metrics.cross_transcript = 0
    metrics.origin_premrna = 0
    metrics.origin_mrna = 0
    metrics.origin_ambiguous = 0

    # The root container:
    data = {}

    LOGGER.info('Parsing BAM file to internal structure...')
    grouped = iCount.mapping.xlsites._processs_bam_file(  # pylint: disable=protected-access
        bam, metrics, mapq_th, strange, annotation=segmentation, gap_th=holesize_th)

    LOGGER.info('Parsing segmentation to internal structure...')
    # pylint: disable=protected-access
    segmentation = iCount.genomes.segment._prepare_annotation(segmentation)

    LOGGER.info('Processing data...')
    for (chrom, strand), by_pos in grouped.items():

        LOGGER.debug("Processing chromosome %s, strand %s", chrom, strand)

        # Sort all genes (and intergenic) by start coordinate.
        segmentation_sorted = sorted(
            segmentation[(chrom, strand)].items(), key=lambda x: x[1]['gene_segment'].start)
        seg_max_index = len(segmentation_sorted) - 1
        start_gene_index, stop_gene_index = 0, seg_max_index

        for xlink_pos, by_bc in sorted(by_pos.items()):
            # pylint: disable=protected-access
            iCount.mapping.xlsites._merge_similar_randomers(by_bc, mismatches)
            # by_bc is modified in place in _merge_similar_randomers

            for _, hits in by_bc.items():
                ss_groups = {}
                for read in hits:
                    ss_groups.setdefault(read[4], []).append(read)

                # Process each group:
                for ss_group in ss_groups.values():

                    # This block extracts just the required genes (& gene_content)
                    # without iterating through all genes/content in chromosome.

                    # Sort reads by length and take the longest one!
                    ss_group = sorted(ss_group, key=lambda x: (-x[2]))
                    stop = ss_group[0][1]
                    start = xlink_pos
                    segmentation_subset = []
                    passed_start = False  # Weather the start_gene_index was aready found.
                    for gene_item in segmentation_sorted[start_gene_index:]:
                        gene_segment = gene_item[1]['gene_segment']

                        if gene_segment.start <= start <= gene_segment.stop:
                            start_gene_index = segmentation_sorted.index(gene_item)
                            passed_start = True

                        if passed_start:
                            segmentation_subset.append(gene_item)

                        if gene_segment.start <= stop <= gene_segment.stop:
                            stop_gene_index = segmentation_sorted.index(gene_item)
                            # Append also one gene before (insert on first position to keep sorted)
                            segmentation_subset.insert(
                                0, segmentation_sorted[max(start_gene_index - 1, 0)])
                            # Append also one gene after:
                            segmentation_subset.append(
                                segmentation_sorted[min(stop_gene_index + 1, seg_max_index)])
                            break
                        # Even if entries repeat, this is still OK, sice
                        # first and last gene neeed to be the ones not including
                        # start/stop!

                    # segmentation_subset is defined. Now process this group:
                    _process_read_group(
                        xlink_pos, chrom, strand, ss_group[0], data, segmentation_subset, metrics,
                        implicit_handling=implicit_handling)

    # Write results (``data``) to file:
    header = ['RNAmap type', 'position', 'all', 'explicit']
    cross_tr_header = ['chrom', 'strand', 'xlink', 'second-start', 'end-position', 'read_len']
    LOGGER.info('Writing output files...')
    with open(out_file, 'wt') as ofile, open(cross_transcript, 'wt') as vfile:
        ofile.write('\t'.join(header) + '\n')
        vfile.write('\t'.join(cross_tr_header) + '\n')
        for rna_map_type, positions in sorted(data.items()):
            if rna_map_type == 'cross_transcript':
                for (chrom, strand, xlink), read_list in positions.items():
                    for (_, end, read_len, _, second_start) in read_list:
                        vfile.write('\t'.join(map(
                            str, [chrom, strand, xlink, second_start, end, read_len])) + '\n')
            else:
                for position, [all_, explic] in sorted(positions.items()):
                    all_, explic = _f2s(all_, dec=4), _f2s(explic, dec=4)
                    ofile.write('\t'.join(map(str, [rna_map_type, position, all_, explic])) + '\n')
    LOGGER.info('RNA-maps output written to: %s', out_file)
    LOGGER.info(
        'Reads spanning multiple transcripts annotation written to %s', cross_transcript)
    LOGGER.info('Done.')
    return metrics
