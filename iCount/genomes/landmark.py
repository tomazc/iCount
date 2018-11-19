""".. Line to protect from pydocstyle D205, D400.

Landmarks
---------

Landmark generation.
"""
from pybedtools import BedTool, create_interval_from_list

from .constants import RNAMAP_TYPES, TYPE_HIERARCHY

LANDMARKS_FILE = 'regions.gtf.gz'


def get_gene_name(seg1, seg2):
    """Determine gene name to use in landmarks file."""
    level1 = TYPE_HIERARCHY.index(seg1[2])
    level2 = TYPE_HIERARCHY.index(seg2[2])

    if level1 < level2:
        return seg1.attrs.get('gene_name', '.')
    elif level2 < level1:
        return seg2.attrs.get('gene_name', '.')
    else:
        raise ValueError('Not possible')


def make_single_type_landmarks(regions, maptype):
    """Create landmarks of single maptype."""
    up_types = RNAMAP_TYPES[maptype]['upstream-type']
    up_limit = RNAMAP_TYPES[maptype].get('upstream-size-limit', 0)
    down_types = RNAMAP_TYPES[maptype]['downstream-type']
    down_limit = RNAMAP_TYPES[maptype].get('downstream-size-limit', 0)

    intervals = []
    for strand in ['+', '-']:
        seg1 = None
        for seg2 in BedTool(regions).filter(lambda x: x.strand == strand):  # pylint: disable=cell-var-from-loop
            if seg1 and seg1.chrom == seg2.chrom:

                # pylint: disable=unsubscriptable-object
                if (strand == '+' and
                        seg1[2] in up_types and
                        seg2[2] in down_types and
                        len(seg1) > up_limit and
                        len(seg2) > down_limit):
                    # BED6 iz zero based and does not include stop nuclotide
                    # Landmark position is the first nucleotide in the downstream feature
                    # + strand: ---up---||---down---
                    #           01234567890123456789
                    # - strand: ---down---||---up---
                    gene_name = get_gene_name(seg1, seg2)
                    name = '{};{}'.format(maptype, gene_name)
                    intervals.append([seg1.chrom, seg1.stop, seg1.stop + 1, name, '.', strand])

                elif (strand == '-' and
                      seg1[2] in down_types and
                      seg2[2] in up_types and
                      len(seg1) > down_limit and
                      len(seg2) > up_limit):
                    gene_name = get_gene_name(seg1, seg2)
                    name = '{};{}'.format(maptype, gene_name)
                    intervals.append([seg1.chrom, seg1.stop - 1, seg1.stop, name, '.', strand])
                # pylint: enable=unsubscriptable-object

            seg1 = seg2

    return intervals


def make_landmarks(regions, landmarks):
    """Create landmarks file from regions file."""
    intervals = []
    for maptype in RNAMAP_TYPES:
        intervals.extend(make_single_type_landmarks(regions, maptype))

    BedTool(create_interval_from_list(list_) for list_ in intervals).sort().saveas(landmarks)
    return landmarks
