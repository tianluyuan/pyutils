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


def colorlines(x, y, ncolors=5, cmapname='copper_r', **kwargs):
    """ Plot a line plot in which the lines change colors
    """
    cmap = plt.get_cmap(cmapname)
    norm = colors.Normalize(vmin=0, vmax=ncolors-1)
    for i in range(ncolors):
        chunksize = len(x)/ncolors
        low = i*chunksize
        # add 1 to keep lines connected
        high = min((i+1)*chunksize+1, len(x))
        plt.plot(x[low:high], y[low:high], color=cmap(norm(i)))
