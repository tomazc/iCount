""".. Line to protect from pydocstyle D205, D400.

Plot heatmap RNA-map
--------------------

Plot heatmap for most covered landmarks.
"""
import csv

import numpy as np
import pandas as pd

from . import rnamap

# pylint: disable=wrong-import-order
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
from matplotlib import pyplot as plt  # pylint: disable=wrong-import-position
# pylint: enable=wrong-import-order


def make_position_to_bin(bins):
    """Create ``position_to_bin`` dictionary."""
    position_to_bin = {}
    for pos in range(bins[0], bins[-1] + 1):
        for bin_start, bin_stop in zip(bins, bins[1:]):
            if bin_start <= pos < bin_stop:
                position_to_bin[pos] = bin_start
                break
        else:
            assert pos == bin_stop
            position_to_bin[pos] = bin_start

    return position_to_bin


def get_top_n_landmarks(fname, top_n, up_limit, down_limit):
    """Get top_n covered landmarks."""
    df, total_cdna = rnamap.parse_results_basic(fname)

    # Create a "sum" column, by which one can sort most covered landmarks.
    df["Sum"] = df.sum(axis=1)
    # Sort and take only top_n covered landmarks. Take index and cast to list.
    top_n_landmarks = df.sort_values(by='Sum', ascending=False).iloc[:top_n].index.tolist()

    return top_n_landmarks, total_cdna


def parse_results(fname, up_limit, down_limit, top_n, nbins=None, binsize=None):
    """Parse RNA-maps results file."""
    # Determine bins.
    if nbins is None and binsize is None:
        raise ValueError('Please define ``binsize`` or ``nbins``.')
    elif nbins and binsize:
        raise ValueError('Either ``binsize`` or ``nbins`` can be defined, but not both.')
    elif nbins:
        nbins = int(nbins)
        bins = list(map(int, np.linspace(up_limit, down_limit, nbins + 1).tolist()))
    elif binsize:
        binsize = int(binsize)
        bins = list(range(up_limit, down_limit, binsize))
        if bins[-1] < down_limit:
            bins.append(down_limit)

    # Create "position_to_bin" dict: a dict that maps every xlink poisition to its bin
    position_to_bin = make_position_to_bin(bins)

    # In first pass, only identify top_n landmarks to be plotted.
    top_n_landmarks, total_cdna = get_top_n_landmarks(fname, top_n, up_limit, down_limit)

    # In second pass, generate landmark-position-score "matrix" for top_n_landmarks
    data = pd.DataFrame(
        data=np.zeros((len(top_n_landmarks), len(bins) - 1)),
        index=top_n_landmarks,
        columns=bins[:-1],
    )
    with open(fname, 'rt') as handle:
        next(handle)  # Skip cdna_info row:

        reader = csv.reader(handle, delimiter='\t')
        header = next(reader)
        for row in reader:
            landmark_name = row[0]
            if landmark_name not in top_n_landmarks:
                continue
            for pos, score in zip(header[1:], row[1:]):
                pos, score = int(pos), int(score)
                if up_limit <= pos <= down_limit:
                    bin_ = position_to_bin[pos]
                    data.at[landmark_name, bin_] = data.at[landmark_name, bin_] + score

    # Perform CPM ormalizaton
    data = data.apply(rnamap.normalize_cpm, args=(total_cdna,))

    return data


def plot_rnaheatmap(fname,
                    outfile=None,
                    up_limit=100,
                    down_limit=100,
                    top_n=100,
                    ax=None,
                    binsize=None,
                    nbins=None,
                    colormap='Greys',
                    ):
    """
    Plot RNA-map heatmap.

    Parameters
    ----------
    fname : str
        RNA-maps result file to plot.
    outfile : str
        Output file.
    up_limit : int
        Upstream plot limit.
    down_limit : int
        Sownstream plot limit.
    top_n : int
        Plot heatmap for top_n best covered landmarks.
    ax : str
        An ``matplotlib.axes.Axes`` instance onto which this plot can
        be drawn. This is useful if you would like to use this function
        to plot this image as a subplot of a more complex figure.
    nbins : int
        Number of bins. Either nbins or binsize can be defined, but not both.
    binsize : int
        Bin size. Either nbins or binsize can be defined, but not both.
    colormap : str
        Colormap to use. Any matplotlib colormap can be used.

    Returns
    -------
    None

    """
    # Make sure limits have correct signs.
    up_limit = -abs(int(up_limit))
    down_limit = abs(int(down_limit))
    top_n = int(top_n)

    # User can provide axes instance (subplot) into which to plot this heatmap.
    # If not given a figure and axes instances are created.
    if ax:
        is_independent = False
        fig = ax.get_figure()
    else:
        is_independent = True
        fig = plt.figure()
        ax = plt.subplot(1, 1, 1)

    data = parse_results(fname, up_limit, down_limit, top_n, binsize=binsize, nbins=nbins)

    heatmap = ax.imshow(
        data,
        aspect='auto',
        cmap=colormap,
        interpolation='none',
        extent=(up_limit, down_limit, top_n, 0),
    )

    ax.set_xlabel('Position')
    ax.set_yticks([])
    ax.set_ylim((top_n, 0))
    ax.set_title('RNA-heatmap')

    fig.colorbar(heatmap, ax=ax, orientation='horizontal')
    if is_independent:
        fig.savefig(outfile)
