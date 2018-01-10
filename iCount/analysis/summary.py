""".. Line to protect from pydocstyle D205, D400.

Cross-link site summary
-----------------------

Report count of cross-link events in each region type.
"""
import logging
import math
import re
import os
import tempfile

from pybedtools import BedTool

import iCount
from iCount.genomes.segment import summary_templates, sort_types_subtypes, TEMPLATE_TYPE, TEMPLATE_SUBTYPE, \
    TEMPLATE_GENE, SUMMARY_TYPE, SUMMARY_SUBTYPE, SUMMARY_GENE

LOGGER = logging.getLogger(__name__)


def summary_reports(annotation, sites, out_dir, templates_dir=None):
    """
    Make summary reports for a cross-link file.

    Parameters
    ----------
    annotation : str
        Annotation file (GTF format). It is recommended to use genome-level segmentation (e.g. regions.gtf.gz), that
        is produced by ``iCount segment`` command.
    sites : str
        Croslinks file (BED6 format). Should be sorted by coordinate.
    out_dir : str
        Output directory.
    templates_dir : str
        Directory containing templates for summary calculation. Made by ``iCount segment`` command. If this argument
        is not provided, summary templates are made on the fly.

    Returns
    -------
    iCount.Metrics
        iCount Metrics object.

    """
    iCount.log_inputs(LOGGER, level=logging.INFO)
    metrics = iCount.Metrics()

    if templates_dir is None:
        templates_dir = tempfile.mkdtemp()
        summary_templates(annotation, templates_dir)

    LOGGER.info('Calculating intersection between cross-link and annotation...')
    # pylint: disable=too-many-function-args,unexpected-keyword-arg
    overlaps = BedTool(sites).intersect(
        BedTool(annotation),
        sorted=True,  # invokes memory efficient algorithm for large files
        s=True,  # only report hits in B that overlap A on the same strand
        wb=True,  # write the original entry in B for each overlap
        nonamecheck=True,  # Do not print warnings about name inconsistency to stdout
    ).saveas()
    # pylint: enable=too-many-function-args,unexpected-keyword-arg
    try:
        overlaps[0]  # will raise Error if overlaps is empty:
    except (IndexError, TypeError):
        raise ValueError('No intersections found. This may be caused by different naming of chromosomes in annotation'
                         'and cross-links file (example: "chr1" vs. "1")')

    type_counter, subtype_counter, gene_counter = {}, {}, {}
    LOGGER.info('Extracting summary data from intersection...')
    for segment in overlaps:
        score = int(segment.score)

        type_ = segment[8]
        type_counter[type_] = type_counter.get(type_, 0) + score

        biotype = re.match(r'.*biotype "(.*?)";', segment[-1])
        biotype = biotype.group(1) if biotype else ''
        biotypes = biotype.split(',')
        for biotype in biotypes:
            sbtyp = iCount.genomes.segment.make_subtype(type_, biotype)
            subtype_counter[sbtyp] = subtype_counter.get(sbtyp, 0) + score / len(biotypes)

        gene_id = re.match(r'.*gene_id "(.*?)";', segment[-1])
        gene_id = gene_id.group(1) if gene_id else None
        gene_counter[gene_id] = gene_counter.get(gene_id, 0) + score

    sum_cdna = 0
    for seg in BedTool(sites):
        sum_cdna += int(seg.score)

    def parse_template(template_file):
        """Parse template file."""
        template = {}
        with open(template_file, 'rt') as ifile:
            for line in ifile:
                line = line.strip().split('\t')
                template[line[0]] = line[1:]
        return template

    LOGGER.info('Writing type report...')
    type_template = parse_template(os.path.join(templates_dir, TEMPLATE_TYPE))
    with open(os.path.join(out_dir, SUMMARY_TYPE), 'wt') as out:
        header = ['Type', 'Length', 'cDNA #', 'cDNA %']
        out.write('\t'.join(header) + '\n')
        for type_, cdna in sorted(type_counter.items(), key=lambda x: sort_types_subtypes(x[0])):
            line = [type_, type_template.get(type_, [-1])[0], math.floor(cdna), cdna / sum_cdna * 100]
            out.write('\t'.join(map(str, line)) + '\n')

    LOGGER.info('Writing subtype report...')
    subtype_template = parse_template(os.path.join(templates_dir, TEMPLATE_SUBTYPE))
    with open(os.path.join(out_dir, SUMMARY_SUBTYPE), 'wt') as out:
        header = ['Subtype', 'Length', 'cDNA #', 'cDNA %']
        out.write('\t'.join(header) + '\n')
        for stype, cdna in sorted(subtype_counter.items(), key=lambda x: sort_types_subtypes(x[0])):
            line = [stype, subtype_template.get(stype, [-1])[0], math.floor(cdna), cdna / sum_cdna * 100]
            out.write('\t'.join(map(str, line)) + '\n')

    LOGGER.info('Writing gene report...')
    gene_template = parse_template(os.path.join(templates_dir, TEMPLATE_GENE))
    with open(os.path.join(out_dir, SUMMARY_GENE), 'wt') as out:
        header = ['Gene name (Gene ID)', 'Length', 'cDNA #', 'cDNA %']
        out.write('\t'.join(header) + '\n')
        for gene_id, cdna in sorted(gene_counter.items()):
            gene_name, length = gene_template.get(gene_id, ['', -1])
            if gene_id == '.':
                gene_name = 'intergenic'
            line = ['{} ({})'.format(gene_name, gene_id), length, math.floor(cdna), cdna / sum_cdna * 100]
            out.write('\t'.join(map(str, line)) + '\n')

    LOGGER.info('Done.')
    return metrics
