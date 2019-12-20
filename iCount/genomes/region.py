""".. Line to protect from pydocstyle D205, D400.

Regions
-------

Make regions file and summary templates.

"""
import itertools
import logging
import math
import os
import re

from pybedtools import BedTool, create_interval_from_list

from .constants import TYPE_HIERARCHY, SUBTYPE_GROUPS

LOGGER = logging.getLogger(__name__)

REGIONS_FILE = 'regions.gtf.gz'
TEMPLATE_TYPE = 'template_type.tsv'
TEMPLATE_SUBTYPE = 'template_subtype.tsv'
TEMPLATE_GENE = 'template_gene.tsv'
SUMMARY_TYPE = 'summary_type.tsv'
SUMMARY_SUBTYPE = 'summary_subtype.tsv'
SUMMARY_GENE = 'summary_gene.tsv'


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
