""".. Line to protect from pydocstyle D205, D400.

Bedgraph conversion
-------------------

Convert from BED6 to bedGraph format.

"""
import logging

import pybedtools

import iCount

LOGGER = logging.getLogger(__name__)


def bed2bedgraph(bed, bedgraph, name='User Track', description='User Supplied Track'):
    """
    Convert from BED6 to bedGraph format.

    For further explanation of ``name`` and ``description`` parameters
    (and also others, that are not yet supported) see:
    https://genome.ucsc.edu/goldenPath/help/customTrack.html#TRACK

    Parameters
    ----------
    bed : str
        Input BED6 file.
    bedgraph : str
        Output bedGraph file.
    name : str
        Track label. Can consist of up to 15 characters.
    description : str
        Track description. Can consist of up to 60 characters.

    Returns
    -------
    None

    """
    iCount.log_inputs(LOGGER, level=logging.INFO)  # pylint: disable=protected-access
    LOGGER.info('Checking input parameters...')
    assert bed.endswith(('.bed', '.bed.gz'))

    header = [
        'track',
        'type=bedGraph',
        'name="{}"'.format(name[:15]),
        'description="{}"'.format(description[:60]),
    ]

    LOGGER.info('Converting %s to %s', bed, bedgraph)
    with open(bedgraph, 'wt') as outfile:
        outfile.write(' '.join(map(str, header)) + '\n')
        for line in pybedtools.BedTool(bed):
            bg_line = [line.chrom, line.start, line.stop, '{}{}'.format(line.strand, line.score)]
            outfile.write('\t'.join(map(str, bg_line)) + '\n')
    LOGGER.info('Done.')
