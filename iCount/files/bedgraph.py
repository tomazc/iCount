""".. Line to protect from pydocstyle D205, D400.

Bedgraph conversion
-------------------

Convert from BED6 to bedGraph format.

"""
import logging
import re

import pybedtools

import iCount

LOGGER = logging.getLogger(__name__)


def bed2bedgraph(bed, bedgraph, name='User Track', description='User Supplied Track', visibility=None, priority=None,
                 color=None, alt_color=None, max_height_pixels=None):
    """
    Convert from BED6 to bedGraph format.

    For further explanation of parameters see:
    https://genome.ucsc.edu/goldenPath/help/customTrack.html#TRACK
    https://genome.ucsc.edu/goldenpath/help/trackDb/trackDbHub.html

    Parameters
    ----------
    bed : str
        Input BED6 file.
    bedgraph : str
        Output bedGraph file.
    name : str
        Track label. Should be shorter than 15 characters.
    description : str
        Track description. Should be shorter than 60 characters.
    visibility: str
        Define the initial display mode of the annotation track. Choose
        among "hide", "dense", "full", "pack" and "squish". Default is
        "dense".
    priority : int
        Defines the track's order relative to other tracks in same group.
    color : str
        Define the main color for the annotation track. The track color
        consists of three comma-separated RGB values from 0-255, e.g
        RRR,GGG,BBB. The default value is 0,0,0 (black).
    alt_color : str
        Allow a color range that varies from color to alt_color.
    max_height_pixels : str
        The limits of vertical viewing space for track, though it is
        configurable by the user. Should be of the format <max:default:min>.

    Returns
    -------
    None

    """
    iCount.log_inputs(LOGGER, level=logging.INFO)  # pylint: disable=protected-access
    LOGGER.info('Checking input parameters...')
    assert bed.endswith(('.bed', '.bed.gz'))
    assert bedgraph.endswith('.bedgraph')

    header = [
        'track',
        'type=bedGraph',
        'name="{}"'.format(name),
        'description="{}"'.format(description),
    ]
    if visibility:
        assert re.fullmatch(r'hide|dense|full|pack|squish', visibility)
        header.append('visibility={}'.format(visibility))
    if priority:
        assert re.fullmatch(r'\d+', str(priority))
        header.append('priority={}'.format(priority))
    if color:
        assert re.fullmatch(r'\d{1,3},\d{1,3},\d{1,3}', color)
        header.append('color={}'.format(color))
    if alt_color:
        assert re.fullmatch(r'\d{1,3},\d{1,3},\d{1,3}', alt_color)
        header.append('altColor={}'.format(alt_color))
    if max_height_pixels:
        assert re.fullmatch(r'\d{1,3}:\d{1,3}:\d{1,3}', max_height_pixels)
        header.append('maxHeightPixels={}'.format(max_height_pixels))

    LOGGER.info('Converting %s to %s', bed, bedgraph)
    with open(bedgraph, 'wt') as outfile:
        outfile.write(' '.join(map(str, header)) + '\n')
        for line in pybedtools.BedTool(bed):
            bg_line = [line.chrom, line.start, line.stop, '{}{}'.format(line.strand, line.score)]
            outfile.write('\t'.join(map(str, bg_line)) + '\n')
    LOGGER.info('Done.')
