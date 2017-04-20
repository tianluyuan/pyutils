import functools
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import colors
import matplotlib.pyplot as plt


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


def restrict_axes(ax, xlim=None, ylim=None):
    """ Given a matplotlib axis *ax*, restricts the axis limits to xlim
    and ylim if the current bounds exceed the limits. Otherwise, leave alone.

    *xlim* and *ylim* should be passed as tuples
    """
    def new_lims(curr_lim, bounds):
        """checks whether the current limit exceeds the bounds and returns
        the appropriate new limits based on bounds or the current
        limit. Reverse order (i.e. left > right) is allowed and accounted for.
        """
        lb, rt = curr_lim
        if curr_lim[0] < curr_lim[1]:
            # normal ordering
            if curr_lim[0] < bounds[0]:
                lb = bounds[0]
            if curr_lim[1] > bounds[1]:
                rt = bounds[1]
        elif curr_lim[0] > curr_lim[1]:
            # reverse ordering
            if curr_lim[0] > bounds[0]:
                lb = bounds[0]
            if curr_lim[1] < bounds[1]:
                rt = bounds[1]

        return lb, rt
    
    if xlim is not None:
        ax.set_xlim(*new_lims(ax.get_xlim(), xlim), auto=None)
    if ylim is not None:
        ax.set_ylim(*new_lims(ax.get_ylim(), ylim), auto=None)
        

