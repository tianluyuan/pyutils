from matplotlib.backends.backend_pdf import PdfPages
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
