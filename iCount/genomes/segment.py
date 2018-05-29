""".. Line to protect from pydocstyle D205, D400.

Segmentation
------------

Parse annotation file into internal iCount structure - segmentation, which is used in almost all further analyses.

Currently, only annotations from ENSEMBl and GENCODE are supported.
http://www.gencodegenes.org/
http://www.ensembl.org

Annotation is composed of three levels (gene level, transcript level and segment level).
Example of annotation (for simplicity, only interval level of transcript1 is shown)::

    Gene level:    |--------------gene1--------------|   |-intergenic-|
                                 |---------gene2--------|

    Transcript l.: |----------transcript1---------|
                           |-------transcript2-------|
                                 |------transcript3-----|

    Segment level: |-CDS-||-intron-||-CDS-||-UTR3-|


Two "versions" of segmentation are produced: transcript-wise-segmentation and genome-wise-segmentation:

    * In transcript-wide segmentation, each transcript is partitioned into intervals. Intervals must span the whole
      transcript, but should not intersect with each other inside transcript. However, higher hierarchy levels:
      transcripts and genes can still intersect each other. As a result, intervals from different genes/transcripts can
      intersect. Intervals in transcript wise segmentation are called segments and the file is called segmentation.

    * In genome-wide segmentation, whole genome is partitioned into intervals. Such intervals must span the whole
      genome, and should also not intersect with each other (neither the ones from different genes/transcripts).
      Intervals in genome wise segmentation are called regions and the file is also called regions.

It is best to present both segmentations and their relation to annotation visualy. Example of annotation::

            ------------------------------------------------------------------------------------->
                          |-----------gene1(G1)-----------|
                          |--------transcript1(A)---------|
                          |-exon--|         |----exon-----|
                                   |------------------gene2(G2)--------------------|
                                   |-----------------transcript2(B)----------------|
                                   |-exon--|        |----exon----|          |-exon-|

Example of (transcript-wise) segmentation. Intron and intergenic intervals are made. Also, exons are converted in
CDS/UTR3/UTR5 or ncRNA::

            ------------------------------------------------------------------------------------->
            |-intergenic-|
                          |-----------gene1(G1)-----------|
                          |--------transcript1(A)---------|
                          |-UTR5--||-intron||-----CDS-----|
                                   |------------------gene2(G2)--------------------|
                                   |-----------------transcript2(B)----------------|
                                   |-UTR5--||intron||-----CDS----||-intron-||-UTR3-|
                                                                                    |-intergenic-|

Example of regions (genome-wise segmentation). Now the annotation is "flat": each nuclotide has one and only one
region. How does one decide which region to keep if there are more overlaping segments? The following hierarchy is
taken into account: CDS > UTR3 > UTR5 > ncRNA > intron > intergenic::

            ------------------------------------------------------------------------------------->
            |-intergenic-||--UTR5-||--UTR5-||-----CDS-----||-CDS-||-intron-||-UTR3-||-intergenic-|

"""
import itertools
import logging
import math
import os
import re
import shutil
import tempfile
from collections import Counter, OrderedDict

from pybedtools import BedTool, create_interval_from_list

import iCount

LOGGER = logging.getLogger(__name__)

REGIONS_FILE = 'regions.gtf.gz'
TEMPLATE_TYPE = 'template_type.tsv'
TEMPLATE_SUBTYPE = 'template_subtype.tsv'
TEMPLATE_GENE = 'template_gene.tsv'
SUMMARY_TYPE = 'summary_type.tsv'
SUMMARY_SUBTYPE = 'summary_subtype.tsv'
SUMMARY_GENE = 'summary_gene.tsv'

TYPE_HIERARCHY = [
    'CDS',
    'UTR3',
    'UTR5',
    'ncRNA',
    'intron',
    'intergenic',
]
SUBTYPE_GROUPS = OrderedDict([
    ('mRNA', [
        'IG_C_gene',
        'IG_D_gene',
        'IG_J_gene',
        'IG_LV_gene',
        'IG_V_gene',
        'TR_C_gene',
        'TR_D_gene',
        'TR_J_gene',
        'TR_V_gene',
        'non_stop_decay',
        'nonsense_mediated_decay',
        'polymorphic_pseudogene',
        'protein_coding',
        'retained_intron',
        'sense_intronic',
        'sense_overlapping',
        '3prime_overlapping_ncRNA',
    ]),
    ('lncRNA', [
        'IG_C_pseudogene',
        'IG_J_pseudogene',
        'IG_V_pseudogene',
        'IG_pseudogene',
        'TR_J_pseudogene',
        'TR_V_pseudogene',
        'TEC',
        'antisense',
        'antisense_RNA',
        'bidirectional_promoter_lncRNA',
        'lincRNA',
        'ncRNA sRNA',
        'non_coding',
        'macro_lncRNA',
        'misc_RNA',
        'misc_RNA_pseudogene',
        'processed_pseudogene',
        'processed_transcript',
        'pseudogene',
        'transcribed_processed_pseudogene',
        'transcribed_unitary_pseudogene',
        'transcribed_unprocessed_pseudogene',
        'translated_processed_pseudogene',
        'unitary_pseudogene',
        'unprocessed_pseudogene',
    ]),
    ('miRNA', ['miRNA', 'miRNA_pseudogene']),
    ('mt_tRNA', ['Mt_tRNA', 'Mt_tRNA_pseudogene']),
    ('mt_rRNA', ['Mt_rRNA']),
    ('rRNA', ['rRNA', 'rRNA_pseudogene']),
    ('sRNA', ['sRNA', 'ribozyme', 'scRNA', 'scRNA_pseudogene', 'vaultRNA']),
    ('snRNA', ['snRNA', 'snRNA_pseudogene']),
    ('snoRNA', ['snoRNA', 'snoRNA_pseudogene', 'scaRNA']),
    ('tRNA', ['tRNA_pseudogene']),
    ('intergenic', ['intergenic']),
])


def construct_borders(seg_filtered):
    """
    Make BED6 file with all possible borders in ``seg_filtered``.

    Image is best explanation. This functions creates "boxes"::

                |-intergenic-|
                |            ||-----------gene1(G1)-----------|
                |            ||--------transcript1(A)---------|
                |            ||-UTR5--||-intron||-----CDS-----|
                |            ||       ||------------------gene2(G2)--------------------|
                |            ||       ||-----------------transcript2(B)----------------|
                |            ||       ||-UTR5--||intron||-----CDS----||-intron-||-UTR3-|
                ----------------------------------------------------------------------->
        boxes:  |            ||       ||       ||      ||     ||     ||        ||      |


    Parameters
    ----------
    seg_filtered : pybedtools.BedTool
        Segmentation as pybedtools.BedTool object. Should exclude entries of type "gene" and
        "transcript".

    Returns
    -------
    str
        Absolute path to BED6 file with borders.

    """
    borders_data = {}
    for seg in BedTool(seg_filtered):
        borders_data.setdefault((seg.chrom, seg.strand), set()).update([seg.start, seg.stop])

    intervals = []
    for (chrom, strand), borders in borders_data.items():
        borders = list(sorted(borders))
        for start, stop in zip(borders, borders[1:]):
            intervals.append(create_interval_from_list([chrom, start, stop, '.', '.', strand]))

    borders_bed = BedTool(interval for interval in intervals).sort().saveas()
    return os.path.abspath(borders_bed.fn)


def simplify_biotype(type_, biotype):
    """Return generalized (broader category) biotype."""
    # First handle 'special' cases:
    if biotype in SUBTYPE_GROUPS['mRNA'] and type_ == 'ncRNA':
        return 'lncRNA'
    if biotype in SUBTYPE_GROUPS['mRNA'] and type_ == 'intron':
        return 'pre-mRNA'

    for group, biotypes in SUBTYPE_GROUPS.items():
        if biotype in biotypes:
            return group

    return biotype


def make_uniq_region(seg, types, biotypes, genes):
    """Make pybedtools.Interval representing unique region."""
    assert len(types) == len(biotypes) == len(genes)

    # In case biotype is '3prime_overlapping_ncRNA', make sure type is UTR3
    for i, biotype in enumerate(biotypes):
        if biotype == '3prime_overlapping_ncRNA':
            types[i] = 'UTR3'

    idxs = None
    # Only consider segments with highest rated type:
    for region_type in TYPE_HIERARCHY:
        if region_type in types:
            idxs = [i for i, typ in enumerate(types) if typ == region_type]
            break
    assert idxs

    # Simplify biotypes and pick unique ones
    biotype_groups = set()
    for i, item in enumerate(biotypes):
        if i in idxs and item is not None:
            # pylint: disable=undefined-loop-variable
            biotype_groups.add(simplify_biotype(region_type, item))
            # pylint: enable=undefined-loop-variable

    # Note that each entry in `genes` is a tuple of form (gene_id, gene_name, gene_size)
    genes = list(set([item for i, item in enumerate(genes) if i in idxs]))
    if len(genes) > 1:
        # In case there are two or more genes, pick the longest one:
        gene_sizes = [gsize for (_, _, gsize) in genes]
        max_index = gene_sizes.index(max(gene_sizes))
        genes = [genes[max_index]]
    assert len(genes) == 1

    attrs = 'gene_id "{}"; gene_name "{}"; biotype "{}";'.format(
        genes[0][0],
        genes[0][1],
        ','.join(sorted(biotype_groups)),
    )
    # pylint: disable=undefined-loop-variable
    return create_interval_from_list(
        [seg.chrom, '.', region_type, seg.start + 1, seg.stop, '.', seg.strand, '.', attrs])
    # pylint: enable=undefined-loop-variable


def merge_regions(nonmerged, out_file):
    """Merge adjacent regions if they have same name (e.g. type, gene and biotypes)."""
    # Sort by chrom, strand, start
    nonmerged_data = sorted([itr for itr in BedTool(nonmerged)], key=lambda x: (x.chrom, x.strand, x.start))
    merged_data = []

    def check_merge(itr):
        """Extract data needed to decide if intervals can be merged."""
        return (itr.chrom, itr.strand, itr[2], itr.attrs.get('biotype'), itr.attrs.get('gene_id'))

    # But merge if the name is the same (type, biotypes and gene are all stored in name):
    for _, group in itertools.groupby(nonmerged_data, key=check_merge):
        # Sort all intervals with same (chrom, strand and name) by start
        ints = sorted(list(group), key=lambda x: x.start)
        # Create "merged" interval: take start of first element in `ints` and stop from the last
        # element, copy other content:
        merged_data.append(create_interval_from_list(ints[0][:4] + [ints[-1][4]] + ints[0][5:]))

    BedTool(itr for itr in merged_data).sort().saveas(out_file)


def get_gene_sizes(segmentation):
    """Calculate gene size for each gene in segmentation."""
    gene_sizes = {'.': 0}  # Fill '.' as this is 'gene_id' for intergenic regions
    for seg in BedTool(segmentation).filter(lambda x: x[2] == 'gene'):
        gene_id = seg.attrs.get('gene_id', None)
        gene_sizes[gene_id] = len(seg)
    return gene_sizes


def make_subtype(type_, biotype):
    """Make subtype from type and biotype."""
    if not biotype:
        return type_
    return '{} {}'.format(type_, biotype)


def sort_types_subtypes(entry):
    """Sort (sub)type by the order defined in TYPE_HIERARCHY and SUBTYPE_GROUPS."""
    def get_index(element, list_):
        """Get index of element in list."""
        if element in list_:
            return [list_.index(element)]
        return [len(list_)]

    entry = entry.strip()
    if ' ' not in entry:  # Assume entries without spaces are just types:
        return [get_index(entry, TYPE_HIERARCHY)]

    # This is subtype
    type_, biotype = entry.split(' ')[:2]
    if biotype == 'pre-mRNA':
        biotype = 'mRNA'
    return [get_index(type_, TYPE_HIERARCHY), get_index(biotype, list(SUBTYPE_GROUPS.keys()))]


def summary_templates(annotation, templates_dir):
    """Make summary templates."""
    type_template, subtype_template, gene_template = {}, {}, {}
    for interval in BedTool(annotation):
        length = len(interval)

        type_ = interval[2]
        type_template[type_] = type_template.get(type_, 0) + length

        biotypes = interval.attrs.get('biotype', '').split(',')
        for biotype in biotypes:
            sbtyp = make_subtype(type_, biotype)
            subtype_template[sbtyp] = subtype_template.get(sbtyp, 0) + length / len(biotypes)

        gene_id = interval.attrs.get('gene_id', '')
        gene_name = interval.attrs.get('gene_name', '')
        current_size = gene_template.get(gene_id, ['', 0])[1]
        gene_template[gene_id] = [gene_name, current_size + length]

    # Write type template
    with open(os.path.join(templates_dir, TEMPLATE_TYPE), 'wt') as outfile:
        for type_, length in sorted(type_template.items(), key=lambda x: sort_types_subtypes(x[0])):
            outfile.write('{}\t{}\n'.format(type_, math.floor(length)))

    # Write subtype template
    with open(os.path.join(templates_dir, TEMPLATE_SUBTYPE), 'wt') as outfile:
        for subtype, length in sorted(subtype_template.items(), key=lambda x: sort_types_subtypes(x[0])):
            outfile.write('{}\t{}\n'.format(subtype, math.floor(length)))

    # Write gene template
    with open(os.path.join(templates_dir, TEMPLATE_GENE), 'wt') as outfile:
        for gene_id, (gene_name, length) in sorted(gene_template.items()):
            line = [gene_id, gene_name, str(math.floor(length))]
            outfile.write('\t'.join(map(str, line)) + '\n')


def make_regions(segmentation, out_dir=None):
    """Make regions file (regions.gtf.gz) and summary templates."""
    if out_dir is None:
        out_dir = os.getcwd()
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    seg_filtered = BedTool(segmentation).filter(lambda x: x[2] not in ['transcript', 'gene']).saveas()

    borders = construct_borders(seg_filtered)
    # pylint: disable=unexpected-keyword-arg, too-many-function-args
    overlaps = BedTool(borders).intersect(
        seg_filtered, sorted=True, s=True, wa=True, wb=True, nonamecheck=True).saveas()
    # pylint: enable=unexpected-keyword-arg, too-many-function-args

    intervals, types, biotypes, genes = [], [], [], []
    gene_sizes = get_gene_sizes(segmentation)
    pseg = overlaps[0]  # Set initial value for "previous segment"
    for seg in overlaps:
        if seg.start != pseg.start or seg.stop != pseg.stop or seg.strand != pseg.strand:
            intervals.append(make_uniq_region(pseg, types, biotypes, genes))
            types, biotypes, genes = [], [], []

        types.append(seg[8])  # In overlaps, this is 3rd column of GTF file

        biotype = re.match(r'.*biotype "(.*?)";', seg[-1])
        biotype = biotype.group(1) if biotype else None
        biotypes.append(biotype)

        gene_id = re.match(r'.*gene_id "(.*?)";', seg[-1])
        gene_id = gene_id.group(1) if gene_id else None
        gene_name = re.match(r'.*gene_name "(.*?)";', seg[-1])
        gene_name = gene_name.group(1) if gene_name else None
        genes.append((gene_id, gene_name, gene_sizes[gene_id]))

        pseg = seg

    intervals.append(make_uniq_region(pseg, types, biotypes, genes))
    nonmerged = BedTool(interval for interval in intervals).saveas()

    # Merge intervals where possible
    merged = os.path.join(out_dir, REGIONS_FILE)
    merge_regions(nonmerged, merged)

    # Finally, make templates
    summary_templates(merged, out_dir)


def _first_two_columns(input_file):
    """Keep just first two columns of file."""
    ofile = tempfile.NamedTemporaryFile(mode='wt', delete=False)
    with open(input_file) as ifile:
        for line in ifile:
            col_1_2 = line.strip().split()[:2]
            ofile.write('\t'.join(col_1_2) + '\n')
    ofile.close()
    return os.path.abspath(ofile.name)


def _a_in_b(first, second):
    """Check if interval a is inside interval b."""
    return first.start >= second.start and first.stop <= second.stop


def _get_biotype(interval):
    """
    Get interval biotype.

    Biotype of interval is equal to transcript biotype value if present,
    else gene biotype if present else value in column 2 (index 1). The
    last option can only happen in some of early ENSEMBL releases.

    Transcript_biotype and gene biotype values are determined from attributes:

        * In ENSEMBL annotations, they are stored in
          ``transcript_biotype`` and ``transcript_type`` attributes.
        * In GENCODE annotations, they are stored in
          ``gene_biotype`` and ``gene_type`` attributes.

    Parameters
    ----------
    interval : pybedtools.Interval
        Interval == line in GTf file.

    Returns
    -------
    string
        Biotype of interval.

    """
    if "transcript_biotype" in interval.attrs:
        return interval.attrs['transcript_biotype']
    elif "transcript_type" in interval.attrs:
        return interval.attrs['transcript_type']
    elif "gene_biotype" in interval.attrs:
        return interval.attrs['gene_biotype']
    elif "gene_type" in interval.attrs:
        return interval.attrs['gene_type']
    else:
        return interval[1]


def _add_biotype_value(interval, biotype):
    """Add biotype value to interval."""
    col8 = interval[8] if interval[8] != '.' else ''
    return create_interval_from_list(interval[:8] + [col8 + ' biotype "{}";'.format(biotype)])


def _add_biotype_attribute(gene_content):
    """
    Add `biotype` attribute to all intervals in gene_content.

    Parameters
    ----------
    gene_content_ : dict
        Intervals in gene separated by transcript id.

    Returns
    -------
    dict
        Same gene_content_ object with added `biotype` attributes.

    """
    gene_content = gene_content.copy()

    # Determine gene biotype:
    gbiotype = _get_biotype(gene_content['gene'])
    # List to keep track of all possible biotypes in gene:
    gene_biotypes = [gbiotype] if gbiotype else []

    for transcript_id, transcript_intervals in gene_content.items():
        if transcript_id == 'gene':
            continue
        first_exon = [i for i in transcript_intervals if i[2] in ['CDS', 'ncRNA']][0]
        biotype = _get_biotype(first_exon)
        gene_biotypes.append(biotype)

        new_intervals = []
        for interval in transcript_intervals:
            new_intervals.append(_add_biotype_value(interval, biotype))
        gene_content[transcript_id] = new_intervals

    # Finally, make also gene biotype: a list of all biotypes in gene,
    # sorted by frequency. Additionally, another sorting is added to sort
    # by alphabet if counts are equal.
    biotype = ', '.join([i[0] for i in sorted(
        sorted(Counter(gene_biotypes).items()), key=lambda x: x[1], reverse=True)])
    gene_content['gene'] = _add_biotype_value(gene_content['gene'], biotype)

    return gene_content


def _check_consistency(intervals):
    """
    Check that intervals in transcript are consistent with theory.

    The test checks two things:

        * interval[i].stop == interval[i + 1].start - this ensures that
          intervals are covering whole transcript and do not overlap
        * that type of given interval is consistent with the type of
          the preceding one


    Parameters
    ----------
    intervals : list
        Intervals representing transcript.

    Returns
    -------
    None
        None.

    Raises
    ------
    AssertionError
        IN case any of the assert statements fails.

    """
    can_follow = {
        '+': {
            'UTR5': ['intron', 'CDS'],
            'CDS': ['intron', 'UTR3'],
            'intron': ['CDS', 'ncRNA', 'UTR3', 'UTR5'],
            'UTR3': ['intron'],
            'ncRNA': ['intron'],
        },
        '-': {
            'UTR3': ['intron', 'CDS'],
            'CDS': ['intron', 'UTR5'],
            'intron': ['CDS', 'ncRNA', 'UTR3', 'UTR5'],
            'UTR5': ['intron'],
            'ncRNA': ['intron'],
        }
    }
    intervals = intervals.copy()
    try:
        index = next(i for i in range(len(intervals)) if intervals[i][2] == 'transcript')
    except StopIteration:
        raise ValueError("No transcript interval in list of intervals.")
    transcript_interval = intervals.pop(index)
    strand = intervals[0].strand

    intervals = sorted(intervals, key=lambda x: x.start)
    assert transcript_interval.start == intervals[0].start
    assert transcript_interval.stop == intervals[-1].stop
    for first, second in zip(intervals, intervals[1:]):
        assert first.stop == second.start
        assert second[2] in can_follow[strand][first[2]]


def _get_non_cds_exons(cdses, exons, intervals):
    """
    Identify (parts of) exons that are not CDS and classify them.

    The intervals in the output list of this function should have
    "UTR3" or "UTR5" in their third field.

    First, all stop codons and CDS intervals are merged where posssible.
    For each stop codon there two options:
        * touching CDS - merge them
        * NOT touching CDS - rename to CDS
        * stop_codon and CDS overlap completely - ignore stop codon

    For *each* exon, many cases are posssible. Here is a "tree" of all
    possible cases for positive strand:

    * CDS inside exon
        * CDS and exon overlap perfectly  #1 - no UTR here
        * CDS and exon do NOT overlap perfectly
            * cds.start != exon.start & cds.stop == exon.stop
                * exon = utr5 + cds  #2
            * cds.start == exon.start & cds.stop != exon.stop
                * exon = cds + utr3  #3
            * cds.start != exon.start & cds.stop != exon.stop
                * exon = utr5 + cds + utr3  #4
    * CDS not inside exon
        * exon completely UTR, just determine wheather UTR3 or UTR5  #5

    Currently, all 5 cases are covered.


    Parameters
    ----------
    cdses : list
        List of all CDS in transcript group.
    exons : list
        List of all exons in transcript group.
    intervals : list
        List of all intervals in transcript group.

    Returns
    -------
    list
        List of UTR intervals in transcript group.

    """
    utrs = []
    int0 = intervals[0]
    strand = int0.strand
    stop_codons = [i for i in intervals if i.fields[2] == 'stop_codon']
    cdses = cdses.copy()

    # Merge stop_codons with cds where posssible:
    replace_cds_indexes = []
    replace_cdses = []
    new_cdses = []
    for stop_codon in stop_codons:
        touching_cds = next((cds for cds in enumerate(cdses) if
                             stop_codon.stop == cds[1].start or cds[1].stop == stop_codon.start), None)
        # in some GTF files stop_codon is also cds, identify such cases:
        complete_overlap = any(cds.start == stop_codon.start and cds.stop == stop_codon.stop for cds in cdses)
        stop_inside_cds = any(_a_in_b(stop_codon, cds) for cds in cdses)
        if complete_overlap or stop_inside_cds:
            continue  # in this case stop_codon is already turned in cds / inside CDS ...
        elif touching_cds:
            cds_index, cds = touching_cds
            replace_cds_indexes.append(cds_index)
            start = min(cds.start, stop_codon.start) + 1
            stop = max(cds.stop, stop_codon.stop)
            replace_cdses.append(create_interval_from_list(cds[:3] + [start, stop] + cds[5:]))
        else:
            new_cdses.append(create_interval_from_list(stop_codon[:2] + ['CDS'] + stop_codon[3:]))
    for index, cds in zip(replace_cds_indexes, replace_cdses):
        cdses[index] = cds
    cdses.extend(new_cdses)

    for exon in exons:
        if not any(_a_in_b(cds, exon) for cds in cdses):
            # no CDS in exon - completely UTR! Just determine UTR3/UTR5:
            if ((exon.strand == '+' and exon.stop >= max([i.stop for i in cdses])) or
                    (exon.strand == '-' and exon.start <= min([i.start for i in cdses]))):
                mode = "UTR3"
            else:
                mode = "UTR5"
            utrs.append(create_interval_from_list(
                exon[:2] + [mode, exon.start + 1, exon.stop, '.', strand, '.', exon[8]]))

        else:
            # CDS in exons! Identify which one:
            cds = next((c for c in cdses if _a_in_b(c, exon)), None)
            assert cds is not None

            if cds.start != exon.start:
                # UTR in the beggining:
                mode = 'UTR5' if exon.strand == '+' else "UTR3"
                utrs.append(create_interval_from_list(
                    exon[:2] + [mode, exon.start + 1, cds.start, '.', strand, '.', exon[8]]))
            if cds.stop != exon.stop:
                # UTR in the end:
                mode = "UTR3" if exon.strand == '+' else "UTR5"
                utrs.append(create_interval_from_list(
                    exon[:2] + [mode, cds.stop + 1, exon.stop, '.', strand, '.', exon[8]]))

    return cdses, utrs


def _filter_col8(interval, keys=None):
    """Filter the content of 9th column (attributes) in a GTF interval."""
    if keys is None:
        keys = ['gene_id', 'gene_name', 'transcript_id', 'transcript_name']
    return ' '.join(['{} "{}";'.format(key, value) for key, value in sorted(interval.attrs.items()) if key in keys])


def _get_introns(exons):
    """
    Calculate the positions of introns and return them as a list of intervals.

    Introns are intervals located between exons. When constructing them,
    one can copy all the data from exons except:

        * start and stop bp for each intron have to be determined
        * score and frame column are set to undefined value: '.'
        * the name of third column is set to 'intron'
        * content referring to exon is removed from last column


    Parameters
    ----------
    exons : list
        List of exons, sorted by coordinates.

    Returns
    -------
    list
        List of pybedtools interval objects representing introns.

    """
    # start of intron is on exon1.stop + 1
    # stop of intron is on exon2.start (since in GTF: e2.start = e2[3] - 1)
    start_stop = [(e1.stop + 1, e2.start) for e1, e2 in zip(exons, exons[1:])]

    # For introns, keep only a subset of key-value pairs from column 8:
    col8 = _filter_col8(exons[0])

    ex1 = exons[0]
    strand = ex1.strand
    return [create_interval_from_list(ex1[:2] + ['intron', i[0], i[1], '.', strand, '.', col8]) for i in start_stop]


def _process_transcript_group(intervals):
    """
    Process a list of intervals in the transcript group.

    Processed intervals should not overlap and should completely span the
    transcript. Their third column should be one of the following:

        * transcript
        * CDS
        * intron
        * UTR3
        * UTR5
        * ncRNA

    Parameters
    ----------
    intervals : list
        List of intervals in transcript to process.

    Returns
    -------
    list
        Modified list of pybedtools interval objects.

    """
    # Container for interval objects
    container = []

    # Get the interval describing whole transcript (make it if not present)
    index = next((i for i in range(len(intervals)) if intervals[i][2] == 'transcript'), None)
    if index is not None:
        container.append(intervals.pop(index))
    else:
        # Manually create "transcript interval":
        int1 = intervals[0]
        col8 = _filter_col8(int1)

        start = min([i.start for i in intervals])
        stop = max([i.stop for i in intervals])
        container.append(create_interval_from_list(int1[:2] + ['transcript', start + 1, stop] + int1[5:8] + [col8]))

    exons = [i for i in intervals if i[2] == 'exon']
    assert exons

    # Sort exones by exom number (reverse if strand == '-'):
    exons = sorted(exons, key=lambda x: int(x.attrs['exon_number']), reverse=exons[0].strand == '-')
    # Confirm that they are really sorted: get their starts:
    starts = [exon.start for exon in exons]
    assert starts == sorted(starts)

    # Gaps between exons are introns. Apend introns to container container:
    container.extend(_get_introns(exons))

    if not {'CDS', 'start_codon', 'stop_codon'} & {i[2] for i in intervals}:
        # If no CDS/stop_codon/start_codon transcript name should be ncRNA.
        container.extend([create_interval_from_list(i[:2] + ['ncRNA'] + i[3:]) for i in intervals])
    else:
        cdses = [i for i in intervals if i[2] == 'CDS']
        # check that all CDSs are within exons:
        for cds in cdses:
            assert any([_a_in_b(cds, exon) for exon in exons])

        # Determine UTR intervals and new_cds (cds joined with stop codons where possible):
        new_cdses, utrs = _get_non_cds_exons(cdses, exons, intervals)
        container.extend(utrs)
        container.extend(new_cdses)

    # check the consistency of container and make a report if check fails:
    try:
        _check_consistency(container)
    except AssertionError as error:
        print('X' * 80)
        for i in sorted(intervals, key=lambda x: x.start):
            print(i)
        print('*' * 80)
        for i in sorted(container, key=lambda x: x.start):
            print(i)
        raise error

    return container


def _complement(gtf, genome_file, strand, type_name='intergenic'):
    """
    Get the complement of intervals in gtf that have strand == `strand`.

    Required structure of genome_file: first column has to be chromosome
    name and the second column has to be chromosome length. Files produced
    with `samtools faidx` (*.fai file extension) respect this rule.

    Possible options for strand param: '+', '-' and '.'.

    Parameters
    ----------
    gtf : str
        Path to GTF file with content.
    genome_file : str
        Path to genome_file (*.fai or similar).
    strand : string
        Strand for which to compute complement.

    Returns
    -------
    str
        Absolute path to GTF file with complement segments.

    """
    assert(strand in ['+', '-', '.'])

    # Filter the content file by strand
    gtf_strand_only = BedTool(gtf).filter(lambda x: x.strand == strand).saveas()

    # Sorted input is required for complement calculation:
    gtf_sorted = gtf_strand_only.sort(faidx=genome_file).saveas()  # pylint: disable=unexpected-keyword-arg

    # This step will fail if gtf_sorted has different set of chromosomes as they are defined in genome_file:
    intergenic_bed = gtf_sorted.complement(g=genome_file).saveas()

    # intergenic_bed is BED3 file. We need to make it GTF. Note the differences in BED and  GTF:
    # https://pythonhosted.org/pybedtools/intervals.html#bed-is-0-based-others-are-1-based
    # This effectively means:
    # gtf_start = bed_start + 1
    # gtf_stop = bed_stop
    if strand == '+':
        col8 = 'ID "interP%5.5d"; gene_id "."; transcript_id ".";'
    elif strand == '-':
        col8 = 'ID "interN%5.5d"; gene_id "."; transcript_id ".";'
    else:
        col8 = 'ID "interB%5.5d"; gene_id "."; transcript_id ".";'
    gtf = BedTool(
        create_interval_from_list([i[0], '.', type_name, str(int(i[1]) + 1), i[2], '.', strand, '.', col8 % n])
        for n, i in enumerate(intergenic_bed)
    ).saveas()

    return os.path.abspath(gtf.fn)


def _get_gene_content(gtf, chromosomes, report_progress=False):
    """
    Give groups of intervals belonging to one gene (as generator).

    Yielded structure is a dictionary that has key-value pairs:

        * 'gene': interval if type gene
        * 'transcript_id#1': intervals corresponding to transcript_id#1
        * 'transcript_id#2': intervals corresponding to transcript_id#2
        ...

    Parameters
    ----------
    gtf : str
        Path to gtf input file.
    chromosomes : list
        List of chromosomes to consider.
    report_progress : bool
        Show progress.

    Returns
    -------
    dict
        All intervals in gene, separated by transcript_id.

    """
    # Lists to keep track of all already processed genes/transcripts:
    gene_ids = []
    transcript_ids = []

    current_transcript = None
    current_gene = None
    gene_content = {}

    def finalize(gene_content):
        """Procedure before returning group of intervals belonging to one gene."""
        if 'gene' not in gene_content:
            # Manually create "gene interval":
            int1 = next(iter(gene_content.values()))[0]
            col8 = _filter_col8(int1)
            start = min([i.start for j in gene_content.values() for i in j])
            stop = max([i.stop for j in gene_content.values() for i in j])
            gene_content['gene'] = create_interval_from_list(int1[:2] + ['gene', start + 1, stop] + int1[5:8] + [col8])
        return gene_content

    length = BedTool(gtf).count()
    progress, j = 0, 0
    for interval in BedTool(gtf):
        j += 1
        if report_progress:
            new_progress = j / length
            progress = iCount._log_progress(new_progress, progress, LOGGER)  # pylint: disable=protected-access

        if interval.chrom in chromosomes:
            # Segments without 'transcript_id' attributes are the ones that
            # define genes. such intervals are not in all releases.
            if interval.attrs['gene_id'] == current_gene:
                if interval.attrs['transcript_id'] == current_transcript:
                    # Same gene, same transcript: just add to container:
                    gene_content[current_transcript].append(interval)
                else:
                    # New transcript - confirm that it is really a new one:
                    current_transcript = interval.attrs['transcript_id']
                    assert current_transcript not in transcript_ids
                    transcript_ids.append(current_transcript)
                    gene_content[current_transcript] = [interval]

            else:  # New gene!
                # First process old content:
                if gene_content:  # To survive the first iteration
                    yield finalize(gene_content)

                # Confirm that it is really new gene!
                current_gene = interval.attrs['gene_id']
                assert current_gene not in gene_ids
                gene_ids.append(current_gene)

                # Make empty container and classify interval
                gene_content = {}
                if interval[2] == 'gene':
                    gene_content['gene'] = interval
                elif 'transcript_id' in interval.attrs:
                    current_transcript = interval.attrs['transcript_id']
                    assert current_transcript not in transcript_ids
                    transcript_ids.append(current_transcript)
                    gene_content[current_transcript] = [interval]
                else:
                    raise Exception("First element in gene content is neither gene or transcript!")

    # for the last iteration:
    yield finalize(gene_content)


def get_segments(annotation, segmentation, fai, report_progress=False):
    """
    Create GTF file with transcript level segmentation.

    Each line in this file should define one of the following elements:

        * gene
        * transcript
        * CDS
        * UTR3
        * UTR5
        * intron
        * ncRNA
        * intergenic

    Name of third field (interval.fields[2]) should correspond to one
    of theese names. Only consider GTF entries of chromosomes given in
    fai file.

    Parameters
    ----------
    annotation : str
        Path to input GTF file.
    segmentation : str
        Path to output GTF file.
    fai : str
        Path to input genome_file (.fai or similar).
    report_progress : bool
        Show progress.

    Returns
    -------
    str
        Absolute path to output GTF file.

    """
    iCount.log_inputs(LOGGER, level=logging.INFO)
    metrics = iCount.Metrics()
    metrics.genes = 0

    # Container for storing intermediate data
    data = []

    LOGGER.debug('Opening genome file: %s', fai)
    fai = _first_two_columns(fai)
    with open(fai) as gfile:
        chromosomes = [line.strip().split()[0] for line in gfile]

    def process_gene(gene_content):
        """
        Process each group of intervals belonging to gene.

        Process each transcript_group in gene_content, add 'biotype'
        attribute to all intervals and include them in `data`.
        """
        assert 'gene' in gene_content

        for id_, transcript_group in gene_content.items():
            if id_ == 'gene':
                continue
            gene_content[id_] = _process_transcript_group(transcript_group)

        # Add biotype attribute to all intervals:
        gene_content = _add_biotype_attribute(gene_content)

        for id_, transcript_group in gene_content.items():
            if id_ == 'gene':
                continue
            data.extend(transcript_group)
        data.append(gene_content['gene'])

    LOGGER.debug('Processing genome annotation from: %s', annotation)
    for gene_content in _get_gene_content(annotation, chromosomes, report_progress):
        process_gene(gene_content)
        LOGGER.debug('Just processed gene: %s', gene_content['gene'].attrs['gene_id'])
        metrics.genes += 1

    # Produce GTF/GFF file from data:
    gtf = BedTool(i.fields for i in data).saveas()

    LOGGER.info('Calculating intergenic intervals...')
    intergenic_pos = _complement(gtf.fn, fai, '+')
    intergenic_neg = _complement(gtf.fn, fai, '-')

    # Join the gtf, intergenic_pos and intergenic_neg in one file:
    file2 = tempfile.NamedTemporaryFile(delete=False)
    for infile in [gtf.fn, intergenic_pos, intergenic_neg]:
        shutil.copyfileobj(open(infile, 'rb'), file2)
    file2.close()

    file3 = BedTool(file2.name).sort().saveas(segmentation)
    LOGGER.info('Segmentation stored in %s', file3.fn)

    LOGGER.info('Making also gene level segmentation...')
    make_regions(segmentation, out_dir=os.path.dirname(os.path.abspath(segmentation)))
    return metrics


def _prepare_segmentation(seg_file, chrom, strand=None):
    """
    Parse segmentation file to hierarchical structure.

    Utility function to transformn segmentation file to the following
    hierarchical structure::

        segmentation = {
            gene_id#1: {
                'gene_segment': gene_segment,
                transcript_id#1: [transcript, exon1, intron1, exon2, ...],
                transcript_id#2: [transcript, exon1, intron1, exon2, ...],
                ...
            },
            gene_id#2: {},
            ...
        }

    Note that intergenic segments have multiple roles: the same intergenic
    segment has the role of gene, transcript and sub-transcript segment. This
    eases the treatment of intergenic intervals in algorithms that use this
    function (rnamaps, xlsites, ...)

    Parameters
    ----------
    seg_file : str
        Path to GTF file, produces by ``get_segments`` function.

    Returns
    -------
    dict
        Segmentation, wrapped in dict with chrom-strand/gene/transcript levels of
        depth.

    """
    segmentation = {}

    for segment in BedTool(seg_file):
        if segment.chrom != chrom:
            continue
        if strand and segment.strand != strand:
            continue

        if segment[2] == 'gene':
            segmentation.setdefault(segment.attrs['gene_id'], {}). \
                setdefault('gene_segment', segment)

        elif segment[2] == 'intergenic':
            # Make artificial_id from chromosome, strand and start:
            fake_gid = 'G_{}_{}_{}'.format(segment.chrom, segment.strand, segment.start)
            fake_tid = 'T_{}_{}_{}'.format(segment.chrom, segment.strand, segment.start)
            segmentation.setdefault(fake_gid, {})['gene_segment'] = segment
            segmentation[fake_gid][fake_tid] = [segment]

        elif segment[2] == 'transcript':
            segmentation.setdefault(segment.attrs['gene_id'], {}). \
                setdefault(segment.attrs['transcript_id'], []). \
                insert(0, segment)  # Ensure that transcript segment is the first one in list.

        else:  # normal segment
            segmentation.setdefault(segment.attrs['gene_id'], {}). \
                setdefault(segment.attrs['transcript_id'], []). \
                append(segment)

    return segmentation
