""".. Line to protect from pydocstyle D205, D400.

Plot combined (distribution & heatmap) RNA-map
----------------------------------------------

Plot combined (distribution & heatmap) RNA-map.
"""
from . import rnaheatmap, rnamap

# pylint: disable=wrong-import-order
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
from matplotlib import pyplot as plt  # pylint: disable=wrong-import-position
# pylint: enable=wrong-import-order


def plot_combined(fname,
                  outfile,
                  up_limit=100,
                  down_limit=100,
                  top_n=100,
                  smoothing=1,
                  nbins=50,
                  binsize=None,
                  colormap='Greys',
                  ):
    """
    Plot combined (distribution & heatmap) RNA-map.

    Parameters
    ----------
    fname : str
        RNA-maps result file to plot.
    outfile : str
        Output file.
    up_limit : int
        Upstream plot limit.
    down_limit : int
        Downstream plot limit.
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

    Returns
    -------
    None

    """
    fig = plt.figure(figsize=(9, 16))
    gridspec = fig.add_gridspec(2, 1, height_ratios=[1, 4])

    distribution = fig.add_subplot(gridspec[0])
    hist = fig.add_subplot(gridspec[1])

    rnamap.plot_rnamap(
        fname,
        outfile=None,
        up_limit=up_limit,
        down_limit=down_limit,
        ax=distribution,
        smoothing=smoothing,
    )

    rnaheatmap.plot_rnaheatmap(
        fname,
        outfile=None,
        up_limit=up_limit,
        down_limit=down_limit,
        top_n=top_n,
        ax=hist,
        binsize=binsize,
        nbins=nbins,
        colormap=colormap,
    )

    hist.set_title('')

    fig.savefig(outfile)
