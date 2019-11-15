""".. Line to protect from pydocstyle D205, D400.

RNA maps
--------

Perform RNA-maps analysis.
"""
import csv
import logging
import os

import pandas as pd
from pybedtools import BedTool

from iCount.genomes.constants import RNAMAP_TYPES
import iCount  # pylint: disable=wrong-import-position

LOGGER = logging.getLogger(__name__)


def get_single_type_landmarks(landmarks, maptype):
    """Get file with landmarks of only certain type."""
    single_type_landmarks = BedTool(landmarks).filter(lambda x: x.name.split(';')[0] == maptype).saveas()
    return single_type_landmarks.fn


def compute_distances(landmarks, sites, maptype):
    """Compute distances between each xlink and it's closest landmark."""
    # pylint: disable=too-many-function-args,unexpected-keyword-arg
    closest = BedTool(sites).closest(  # pylint: disable=assignment-from-no-return
        landmarks,
        s=True,
        t='first',
        D='a',
        nonamecheck=True,
    )
    # pylint: enable=too-many-function-args,unexpected-keyword-arg

    distances = {}
    total_cdna = 0
    for mseg in closest:
        # "mseg" means merged segment, since it consists of 3 parts:
        # sites segment (BED6) + landmark segment (BED6) + distance
        if len(mseg.fields) != 13:
            LOGGER.warning('Segment length shoud be 13, not %s. Segment fields: %s', str(len(mseg)), str(mseg.fields))

        score = int(mseg[4])
        chrom = mseg[6]
        pos = mseg[7]
        gene_name = mseg[9]
        strand = mseg[11]
        distance = -int(mseg[-1])
        total_cdna += score

        mapdata = RNAMAP_TYPES[maptype]
        if not -mapdata['upstream-size-limit'] <= distance <= mapdata['downstream-size-limit']:
            continue
        if chrom == '.' or strand not in ['+', '-'] or '-' in pos:
            continue

        # loc = Landmark "ID" = landmark exact coordinates
        loc = '{}__{}__{}__{}'.format(chrom, strand, pos, gene_name)

        distances.setdefault(loc, {})[distance] = distances.get(loc, {}).get(distance, 0) + score

    return distances, total_cdna


def make_results_raw_file(distances, fname, total_cdna, maptype):
    """Write distances data to file."""
    up_limit = -RNAMAP_TYPES[maptype]['upstream-size-limit']
    down_limit = RNAMAP_TYPES[maptype]['downstream-size-limit']
    header = list(range(up_limit, down_limit + 1))

    with open(fname, 'wt') as handle:
        outfile = csv.writer(handle, delimiter='\t')
        outfile.writerow(['total_cdna:{}'.format(total_cdna)])

        outfile.writerow(['.'] + header)
        for loc, positions in sorted(distances.items()):
            outfile.writerow([loc] + [positions.get(pos, 0) for pos in header])


def make_results_summarised_file(outdir, fname):
    """Write "plot data" to file."""
    data = {}
    up_limit = -max([mtyp['upstream-size-limit'] for mtyp in RNAMAP_TYPES.values()])
    down_limit = max([mtyp['downstream-size-limit'] for mtyp in RNAMAP_TYPES.values()])
    header = list(range(up_limit, down_limit + 1))

    for basename in [file_ for file_ in os.listdir(outdir) if file_.endswith('.tsv')]:
        results_file = os.path.join(outdir, basename)
        maptype = iCount.plotting.rnamap.guess_maptype(results_file)
        plot_data, _ = iCount.plotting.rnamap.parse_results(results_file)

        # Make date of the same size, impute with 0 score on locations with no data.
        data[maptype] = [plot_data.get(pos, 0) for pos in header]
        assert len(header) == len(data[maptype])

    dframe = pd.DataFrame.from_dict(data, orient='index', columns=header)
    dframe.to_csv(fname, sep='\t')


def run(sites,
        landmarks,
        outdir=None,
        plot_type='combined',
        top_n=100,
        smoothing=1,
        nbins=50,
        binsize=None,
        colormap='Greys',
        imgfmt='png',
        ):
    """
    Compute distribution of cross-links relative to genomic landmarks.

    Parameters
    ----------
    sites : str
        Croslinks file (BED6 format). Should be sorted by coordinate.
    landmarks : str
        Landmark file (landmarks.bed.gz) that is produced by ``iCount segment``.
    outdir : str
        Output directory.
    plot_type : str
        What kind of plot to make. Choices are distribution, heatmaps and combined.
    top_n : int
        Plot heatmap for top_n best covered landmarks.
    smoothing : int
        Smoothing half-window. Average smoothing is used.
    nbins : int
        Number of bins. Either nbins or binsize can be defined, but not both.
    binsize : int
        Bin size. Either nbins or binsize can be defined, but not both.
    colormap : str
        Colormap to use. Any matplotlib colormap can be used.
    imgfmt : str
        Output image format.

    Returns
    -------
    None

    """
    iCount.logger.log_inputs(LOGGER)

    assert plot_type in ['distribution', 'heatmap', 'combined']

    sites_name = iCount.files.remove_extension(sites, ['.bed', '.bed.gz'])
    if outdir is None:
        outdir = os.path.join(os.path.abspath(os.getcwd()), 'rnamaps_{}'.format(sites_name))
        LOGGER.info('Output directory not given, creating one at %s', outdir)
    os.makedirs(outdir, exist_ok=True)

    for maptype, mapdata in RNAMAP_TYPES.items():
        LOGGER.info('Creating landmarks for %s', maptype)
        landmarks_single_type = get_single_type_landmarks(landmarks, maptype)
        if not BedTool(landmarks_single_type).head(n=1, as_string=True):
            LOGGER.info('No landmarks for %s', maptype)
            continue

        LOGGER.info('Processing data for %s', maptype)
        distances, total_cdna = compute_distances(landmarks_single_type, sites, maptype)
        if not distances:
            LOGGER.warning('No distances for %s', maptype)
            continue

        LOGGER.info('Writing results to file for %s', maptype)
        # File with full set of results:
        results_raw_file = os.path.join(outdir, '{}_{}.tsv'.format(sites_name, maptype))
        make_results_raw_file(distances, results_raw_file, total_cdna, maptype)

        kwargs = {
            'fname': results_raw_file,
            'outfile': os.path.join(outdir, '{}_{}.{}'.format(sites_name, maptype, imgfmt)),
            'up_limit': mapdata.get('upstream-plot-width', mapdata['upstream-size-limit']),
            'down_limit': mapdata.get('downstream-plot-width', mapdata['downstream-size-limit']),
            'top_n': top_n,
            'smoothing': smoothing,
            'nbins': nbins,
            'binsize': binsize,
            'colormap': colormap,
        }
        if plot_type == 'distribution':
            kwargs['fnames'] = kwargs['fname']
            for key in ['top_n', 'nbins', 'binsize', 'colormap', 'fname']:
                del kwargs[key]
            iCount.plotting.rnamap.plot_rnamap(**kwargs)  # pylint: disable=unexpected-keyword-arg
        elif plot_type == 'heatmap':
            del kwargs['smoothing']
            iCount.plotting.rnaheatmap.plot_rnaheatmap(**kwargs)  # pylint: disable=unexpected-keyword-arg
        elif plot_type == 'combined':
            iCount.plotting.rnacombined.plot_combined(**kwargs)

    # Single file with only RNA-maps distibution plot data.
    results_summarised_file = os.path.join(outdir, '{}_plot_data.tsv'.format(sites_name))
    make_results_summarised_file(outdir, results_summarised_file)

    LOGGER.info('Done.')
