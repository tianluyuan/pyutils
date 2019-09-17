import functools
import numpy as np
from scipy import optimize
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import colors
import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
from collections import namedtuple


CoordSys = namedtuple('CoordSys', 'det eq gal')
COORD = CoordSys('Detector', 'Equatorial', 'Galactic')


def multipage(filename, figs=None, dpi=200):
    """ http://stackoverflow.com/questions/26368876/saving-all-open-matplotlib-figures-in-one-file-at-once
    """
    pp = PdfPages(filename)
    if figs is None:
        figs = [plt.figure(n) for n in plt.get_fignums()]
    for fig in figs:
        fig.savefig(pp, format='pdf')
    pp.close()


def colorlines(x, y, ncolors=5, cmapname='viridis_r', **kwargs):
    """Plot a line plot in which the lines change colors as the data is
    stepped through.

    *ncolors* specifies the number of different colors to use
    """
    cmap = plt.get_cmap(cmapname)
    norm = colors.Normalize(vmin=0, vmax=ncolors-1)
    for i in range(ncolors):
        chunksize = len(x)/ncolors
        low = i*chunksize
        # add 1 to keep lines connected
        high = min((i+1)*chunksize+1, len(x))
        plt.plot(x[low:high], y[low:high], color=cmap(norm(i)), **kwargs)


def pdf(func):
    """ decorator to save all plots generated by func to pdf
    """
    @functools.wraps(func)
    def pdfwrapper(*args, **kwargs):
        if kwargs.has_key('pdffile') and kwargs['pdffile'] is not None:
            plt.close('all')

            ret = func(*args, **kwargs)
            multipage(kwargs['pdffile'])
            return ret
        else:
            return func(*args, **kwargs)

    return pdfwrapper


def new_lims(curr_lim, bounds):
    """checks whether the current limit exceeds the bounds and returns
    the appropriate new limits based on bounds or the current
    limit. Reverse order (i.e. left > right) is allowed and accounted for.
    """
    lb, rt = curr_lim
    if lb <= rt:
        # normal ordering
        combined = sorted(curr_lim+bounds)
        # no overlap
        if tuple(combined[:2]) == (lb, rt) or tuple(combined[2:]) == (lb, rt):
            return bounds[0], bounds[1]

        return combined[1], combined[2]
    else:
        # reverse ordering
        combined = sorted(curr_lim+bounds, reverse=True)
        # no overlap
        if tuple(combined[:2]) == (lb, rt) or tuple(combined[2:]) == (lb, rt):
            return bounds[0], bounds[1]

        return combined[1], combined[2]
    

def restrict_axes(ax, xlim=None, ylim=None):
    """Given a matplotlib axis *ax*, restricts the axis limits to xlim and
    ylim if they exceed the bounds (xlim, ylim). If the axis limit
    does not overlap with (xlim, ylim), the new limits are set to
    (xlim, ylim). Otherwise limits are kept as is.

    *xlim* and *ylim* are the restricted ranges and should be passed as tuples
    """
    if xlim is not None:
        ax.set_xlim(*new_lims(ax.get_xlim(), xlim), auto=None)
    if ylim is not None:
        ax.set_ylim(*new_lims(ax.get_ylim(), ylim), auto=None)
        

def contour_levels(x, y, cls=(0.95, 0.68), bins=None):
    """given 2D datapoints, return values of the pdf corresponding to the
    passed confidence levels
    """
    if bins is None:
        bins = int(np.sqrt(len(x)))
    # Make a 2d normed histogram
    H,xedges,yedges=np.histogram2d(x,y,bins=bins,normed=True)

    norm=H.sum() # Find the norm of the sum

    # Take histogram bin membership as proportional to Likelihood
    # This is true when data comes from a Markovian process
    def objective(limit, target):
        w = np.where(H>limit)
        count = H[w]
        return count.sum() - target

    levels = [optimize.bisect(objective, H.min(), H.max(), args=(cl*norm,))
              for cl in cls]
    levels.append(H.max())
    return levels


def hp_ticklabels(coord, zoom=False, lonra=None, latra=None, rot=None):
    """ labels coordinates on a healpy map

    zoom: indicates zoomed-in cartview
    lonra: longitude range of zoomed-in map
    latra: latitude range of zoom-in map
    rot: center of zoomed in map
    lcoord: label of coordinate system
    """
    import healpy as hp
    # coordinate labels
    ax = plt.gca()
    if zoom:
        # location of other, fixed coordinate
        lon_offset = rot[0]+lonra[0]
        lat_offset = rot[1]+latra[0]
        # lonlat coordinates for labels
        lons = np.arange(np.round(lon_offset),
                         lon_offset+lonra[1]-lonra[0], 2)
        lats = np.arange(np.round(lat_offset),
                         lat_offset+latra[1]-latra[0], 2)
    else:
        lon_offset = -180
        lat_offset = 0

        # lonlat coordinates for labels
        lons = np.arange(-150, 181, 30)
        lats = np.arange(-90, 91, 30)

    # actual text at those coordinates
    if coord == COORD.det:
        llats = 90-lats
    else:
        llats = lats

    # white outline around text
    pe = [path_effects.Stroke(linewidth=1.5, foreground='white'),
          path_effects.Normal()]
    for _ in zip(lats, llats):
        hp.projtext(lon_offset, _[0], "{:.0f}$^\circ$".format(_[1]),
                    lonlat=True, path_effects=pe)
    if zoom:
        for _ in lons:
            hp.projtext(_, lat_offset,
                        "{:.0f}$^\circ$".format(_), lonlat=True,
                        path_effects=pe)
    else:
        ax.annotate(r"$\bf{-180^\circ}$", xy=(1.7, 0.625), size="medium")
        ax.annotate(r"$\bf{180^\circ}$", xy=(-1.95, 0.625), size="medium")
    ax.annotate(coord, xy=(0.8, -0.05),
                size="medium", xycoords="axes fraction")

