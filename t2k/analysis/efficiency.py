from tianlu.rootutils import (glob_to_chain,
                              split_varstring,
                              get_hist)
from tianlu.analysis import utils as tautils
from collections import namedtuple
import ROOT as rt


__author__ = 'Tianlu Yuan'
__email__ = 'tianlu.yuan [at] colorado.edu'


# namedtuple classes for setting up numer and denom and also vars
Options = namedtuple('Options', 'name cut tree_key')
Var = namedtuple('Var', 'var binning title')

class Efficiency:
    def __init__(self, filepath, num_tup, denom_tup, var_tup):
        """ This class helps create a TEfficiency object based on the passed
        arguments.

        filepath: wildcard path or list of paths to root files
        num/denom_tup: namedtuples of Options() containing necessary name, cutstring,
        and tree to use for creating the numerator and denominator histos
        var_tup: namedtuple of Var() that contains variable to plot the efficiency
        as a function of, binning, and title
        """
        num_chain = rt.TChain(num_tup.tree_key)
        glob_to_chain(num_chain, filepath)

        if denom_tup.tree_key == num_tup.tree_key:
            denom_chain = num_chain
        else:
            denom_chain = rt.TChain(denom_tup.tree_key)
            glob_to_chain(denom_chain, filepath)

        var = var_tup.var
        binning = var_tup.binning

        self.h_numer = get_hist('h_numer',
                                num_tup.name,
                                binning)
        self.h_numer.Sumw2()
        self.h_denom = get_hist('h_denom',
                                num_tup.name,
                                binning)
        self.h_denom.Sumw2()

        num_chain.Draw(var+">>h_numer", num_tup.cut)
        denom_chain.Draw(var+">>h_denom", denom_tup.cut)

        self.efficiency = rt.TEfficiency(self.h_numer, self.h_denom)
        rt.SetOwnership(self.efficiency, 0)
        self.efficiency.SetLineWidth(1)

        # this should but doesn't set the axes in ROOT 5.34/09
        # frac = '{0}/{1}'.format(num_tup.name, denom_tup.name)
        frac = 'Fraction'
        title = '{0};{1};{2}'.format(num_tup.name,
                                     var.replace(':', ';'),
                                     frac)
        self.efficiency.SetTitle(title)
        # workaround
        self.set_axes(var, frac)


    def set_axes(self, var, frac):
        """ Alternative method to setting axes titles for ROOT 5.34/09
        """
        self.efficiency.Draw()
        rt.gPad.Update()

        # use regex to match varstrings for two dimensions'
        split_vars = split_varstring(var)
        dim = len(split_vars)

        if dim == 1:
            tautils.set_axes_titles(var,
                                    self.efficiency.GetPaintedGraph(),
                                    frac)
            self.efficiency.GetPaintedGraph().GetYaxis().SetRangeUser(0,1)
        elif dim == 2:
            tautils.set_axes_titles(var,
                                    self.efficiency.GetPaintedHistogram(),
                                    frac)
            self.efficiency.GetPaintedHistogram().GetZaxis().SetRangeUser(0,1)


    def get_legend_text(self):
        # create legend with overall ratios
        n_numer = self.h_numer.Integral(0, self.h_numer.GetNbinsX()+1)
        n_denom = self.h_denom.Integral(0, self.h_denom.GetNbinsX()+1)
        legend_text = (self.efficiency.GetTitle()+' (overall: %.4f)'
                       % (float(n_numer)/n_denom))

        print legend_text
        return legend_text


    def set_color(self, idx):
        self.efficiency.SetMarkerColor(idx)
        self.efficiency.SetLineColor(idx)


    def set_style(self, idx):
        self.efficiency.SetLineStyle(idx)


    def draw(self, opt=''):
        self.efficiency.Draw(opt)
        rt.gPad.Update()


    def draw_numerator(self, opt=''):
        self.h_numer.Draw(opt)


    def draw_denominator(self, opt=''):
        self.h_denom.Draw(opt)
