""" Useful ROOT utilities
"""
import os
import sys
import glob
import re
import math
import functools
from collections import namedtuple
from array import array
import ROOT as rt
from common import utils


__author__ = 'Tianlu Yuan'
__email__ = 'tianlu.yuan [at] colorado.edu'


GOODCOLORS = [1, 2, 4, 797, 904, 8, 9, 882, 630, 598, 792]

# Create a namedtuple, MultiCanvas object.  ROOT's handling of
# manually drawn pads makes it difficult to retrieve the pad after
# it's drawn on the canvas.  The MultiCanvas object will make things
# more pythonic
class MultiCanvas(object):
    def __init__(self, canv, pads):
        self.canvas = canv
        self.pads = pads


    def Close(self):
        self.canvas.Close()


    def Clone(self):
        return self.canvas.Clone()


    def Update(self):
        for pad in self.pads:
            pad.Update()
        self.canvas.Update()



class Efficiency:
    # namedtuple classes for setting up numer and denom and also vars
    Options = namedtuple('Options', 'name cut tree_key')
    Var = namedtuple('Var', 'var binning')

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
        etitle = '{0};{1};{2}'.format(num_tup.name,
                                      var.replace(':', ';'),
                                      frac)
        self.efficiency.SetTitle(etitle)


    def get_legend_text(self):
        # create legend with overall ratios
        n_numer = self.h_numer.Integral(0, self.h_numer.GetNbinsX()+1)
        n_denom = self.h_denom.Integral(0, self.h_denom.GetNbinsX()+1)
        legend_text = (self.efficiency.GetTitle()+' (overall: %.4f)'
                       % (float(n_numer)/n_denom))

        print(legend_text)
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


# convert ROOT bool to python bool (from robj)
# deprecated in python 2.6 with ROOT v5.34
def convert_bool(aROOTBool):
    return bool(int(aROOTBool.encode('hex')))


def set_stack_range(stack, lim='max', opt=''):
    """Sets the max of THStack with or without the option "noStack"

    lim: determines whether to set the max, min, or both based on the
    passed THStack stack. There is also an option 'ratio', which sets
    a symmetric range about 1.0 based on the maximum.
    """
    stackMax = stack.GetMaximum(opt)
    stackMin = stack.GetMinimum(opt)
    stackRange = stackMax - stackMin

    if lim == 'max':
        stack.SetMaximum(stackMax + 0.1 * stackRange)
        return
    elif lim == 'min':
        stack.SetMinimum(stackMin - 0.1 * stackRange)
        return
    elif lim == 'both':
        stack.SetMaximum(stackMax + 0.1 * stackRange)
        stack.SetMinimum(stackMin - 0.1 * stackRange)
        return
    elif lim == 'ratio':
        stack.SetMaximum(stackMax + 0.1 * stackRange)
        stack.SetMinimum(2 - stackMax - 0.1 * stackRange)
        return
    else:
        sys.exit('Error setting stack free axis limits.  Invalid option.')


def get_legend(x1=0.40, x2=0.90, y1=0.75, y2=0.95):
    leg = rt.TLegend()
    leg.SetX1NDC(x1)
    leg.SetX2NDC(x2)
    leg.SetY1NDC(y1)
    leg.SetY2NDC(y2)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    return leg


def loop_bins(func):
    ''' Decorator to loop over bins in a histogram and run some function on 
    each bin.
    '''
    @functools.wraps(func)
    def loop_func(h_temp):
        for i in range(h_temp.GetNbinsX()+2):
            for j in range(h_temp.GetNbinsY()+2):
                for k in range(h_temp.GetNbinsZ()+2):
                    func(h_temp, i, j, k)
    return loop_func


@loop_bins
def unit_hist(h_temp, i, j, k):
    h_temp.SetBinContent(i,j,k, 1.0)
    h_temp.SetBinError(i,j,k, 0.0)


@loop_bins
def reset_errors(h_temp, i, j, k):
    h_temp.SetBinError(i,j,k, 0)


@loop_bins
def reset_poisson_errors(h_temp, i, j, k):
    h_temp.SetBinError(i,j,k, math.sqrt(h_temp.GetBinContent(i,j,k)))


def project_hist(hist):
    """ Hist the histogram's bin contents
    """
    nbins = int(math.sqrt(hist.GetSize()))
    minx = hist.GetMinimum() - 0.1 * abs(hist.GetMinimum())
    maxx = hist.GetMaximum() + 0.1 * abs(hist.GetMaximum())

    hproj = rt.TH1F('hproj', 'Projection', nbins, minx, maxx)
    hproj.GetXaxis().SetTitle(get_primary_axis(hist).GetTitle())
    hproj.GetYaxis().SetTitle('N')
    [hproj.Fill(hist.GetBinContent(i)) for i in range(hist.GetSize())
     if not hist.IsBinOverflow(i) and not hist.IsBinUnderflow(i)]
    return hproj


def scale_multilabel(h1, h2, p1, p2):
    """Histgram axes have sizes as a fraction of the pad width (default is
    0.04). For multicanvases this leads different label sizes. This
    function scales the histograms labels to the larger size.
    """
    ratio = float(p1.GetWh()*p1.GetAbsHNDC())/(p2.GetWh()*p2.GetAbsHNDC())
    h_small, scale = (h1, 1/ratio) if ratio < 1 else (h2, ratio)

    orig_scale = h_small.GetXaxis().GetLabelSize()
    orig_offset = h_small.GetXaxis().GetTitleOffset()

    h_small.GetXaxis().SetLabelSize(orig_scale*scale)
    h_small.GetYaxis().SetLabelSize(orig_scale*scale)
    h_small.GetXaxis().SetTitleSize(orig_scale*scale)
    h_small.GetYaxis().SetTitleSize(orig_scale*scale)
    # h_small.GetXaxis().SetTitleOffset(orig_offset/scale)
    h_small.GetYaxis().SetTitleOffset(orig_offset/scale)


def get_primary_axis(hist):
    """Returns the highest dimension axis, e.g. if 1D histo this will
    return the y-axis.
    """
    dim = hist.GetDimension()
    if dim == 1:
        return hist.GetYaxis()
    elif dim == 2:
        return hist.GetZaxis()


def label_primary_axis(hist, label):
    """Labels the highest dimension axis, i.e. if 1D histo this will
    label the y-axis.
    """
    paxis = get_primary_axis(hist)
    if paxis:
        paxis.SetTitle(label)
    else:
        hist.SetTitle(label)


def residuals(h1, h2, normalized=True, chi2opt='UU'):
    """Returns histogram of the residuals between h1 and h2. If
    normalized is true the residuals are calculated via the
    TH1::Chi2Test
    """
    hres = h1.Clone()
    # calculate the variance summed across all bins. If this is zero,
    # then cannot get normalized residuals
    totvar = sum([hres.GetBinError(i)**2 for i in range(hres.GetSize())])

    if not normalized or totvar==0:
        hres.Add(h2, -1)
        label_primary_axis(hres, 'Residuals')
    else:
        res = array('d', [0.0]*h1.GetSize())
        # Chi2Test for multi dimensional histogram gives incorrect
        # residuals, flatten the histogram first and get residuals
        # using flattened version, when drawn it will retain the
        # original number of dimensions. Include overflows.
        h1_flat = flatten_hist(h1, True)
        h2_flat = flatten_hist(h2, True)
        h1_flat.Chi2Test(h2_flat, chi2opt, res)
        [hres.SetBinContent(idx, content) for idx, content in enumerate(res)]
        reset_errors(hres)
        label_primary_axis(hres, 'Normalized Residuals')

    hres.SetFillColor(0)
    hres.SetLineColor(rt.kBlack)
    hres.SetLineStyle(1)
    return hres


def draw_multi_canvas(h_above, h_below, opt_above='', opt_below=''):
    """ Splits main_pad into two partitions (e.g. for residual plots).
    Draws THObject h_above in upper pad and h_below in lower pad.

    returns TCanvas containing both histos

    Taken from residuals.cxx example on root forum
    """
    multi_canvas = get_multi_canvas()

    rt.gStyle.SetOptStat(0)

    pad1 = multi_canvas.pads[0]
    pad2 = multi_canvas.pads[1]
    scale_multilabel(h_above, h_below, pad1, pad2)

    pad1.cd()
    h_above.Draw(opt_above)
    pad2.cd()
    h_below.Draw(opt_below)

    canv_clone = multi_canvas.Clone()
    multi_canvas.Close()
    return canv_clone


def get_multi_canvas():
    """ Splits main_pad into two partitions (e.g. for residual plots).

    returns a MultiCanvas tuple of TCanvas with the two subdivided pads

    Taken from residuals.cxx example on root forum
    """
    canv = rt.TCanvas('canv')
    canv.cd()

    pad1 = rt.TPad('pad1', 'pad1', 0, 0.33, 1, 1)
    pad2 = rt.TPad('pad2', 'pad2', 0, 0, 1, 0.33)

    # some styling
    pad1.SetBottomMargin(0.00001)
    pad1.SetBorderMode(0)
    pad2.SetTopMargin(0.00001)
    pad2.SetBottomMargin(0.2)
    pad2.SetBorderMode(0)

    # draw pad onto current canvas
    pad1.Draw()
    pad2.Draw()

    # create python container for multi-canvas
    multi_canvas = MultiCanvas(canv, [pad1, pad2])

    return multi_canvas




def draw_inset(main_canv, inset_canv):
    """Draws an inset canvas within the main canvas. Both main_canv and
    inset_canv must be TCanvases
    """
    multi_inset = get_inset_canvas()
    multi_inset.pads[0].cd()
    main_canv.DrawClonePad()
    multi_inset.pads[1].cd()
    inset_canv.DrawClonePad()

    input('Enter to continue...')
    inset_canv_clone = multi_inset.canvas.Clone()
    multi_inset.Close()
    return inset_canv_clone


def get_inset_canvas():
    """ Splits main_pad into two with smaller one inset.

    returns a MultiCanvas tuple of TCanvas with the two subdivided pads

    From https://root.cern.ch/phpBB3/viewtopic.php?f=3&t=11578
    """
    canv = rt.TCanvas('canv_inset')
    canv.cd()

    mainpad = rt.TPad('mainpad', 'mainpad', 0, 0, 1, 1)
    subpad = rt.TPad('subpad', 'subpad', 0.5, 0.4, 0.9, 0.8)

    # draw pad onto current canvas
    mainpad.Draw()
    subpad.Draw()

    # create python container for multi-canvas
    multi_canvas = MultiCanvas(canv, [mainpad, subpad])

    return multi_canvas


def vector(cptype, python_list):
    """ Converts a python_list to a std::vector.

    cptype: the c++ type of the vector.  Should be either a string of the name
    of the type or of the form ROOT.TObject

    python_list: a list of the elements that will be pushed onto the vector.
    This list is not type-checked

    returns: ROOT.std.vector object
    """
    v = rt.std.vector(cptype)()

    for ele in python_list:
        v.push_back(ele)

    return v


def get_fracerrs(hist):
    """ hist is a ROOT.TH object.  Currently implemented for up to 3D hist.

    returns: TH object of fractional errors, cloned from hist
    """
    h_ferror = hist.Clone()
    h_ferror.GetYaxis().SetTitle('Fractional Error')

    @loop_bins
    def errors_single_bin(hist, i, j, k):
        this_bin = hist.GetBin(i, j, k)
        bin_content = hist.GetBinContent(this_bin)
        bin_error = hist.GetBinError(this_bin)

        if bin_content != 0 and this_bin != 0:
            frac_error = bin_error/bin_content
        else:
            frac_error = 0.

        h_ferror.SetBinContent(this_bin, frac_error)
        h_ferror.SetBinError(this_bin, 0)

    errors_single_bin(hist)

    return h_ferror


def get_fracbias(meas, true):
    """meas, and true are ROOT.TH objects with identical
    binning. Currently implemented for up to 3D hist.

    returns: TH object of fractional bias (meas-true)/true, cloned from hist

    """
    bias = meas.Clone()
    bias.GetYaxis().SetTitle('Fractional Bias')

    bias.Add(true, -1)
    bias.Divide(true)
    reset_errors(bias)

    return bias


def get_branch_ctype(branch):
    """ Returns the branch variable C type.
    """
    if isinstance(branch, rt.TBranchElement):
        # TBranchElement is branch for object types
        # it contains extra members to get the type
        return branch.GetTypeName()
    elif isinstance(branch, rt.TBranch):
        branch_split = branch.GetTitle().split('/')
        if len(branch_split) > 1:
            return branch_split[1]
        else:
            raise RuntimeError('Branch "{0}" has no ctype!'.
                               format(branch.GetName()))
    else:
        raise TypeError('Branch "{0}" is not of type "TBranch".'.
                        format(branch.GetName()))



def get_branch_pytype(branch):
    """ Returns the branch type for use in python arrays.
    Converts from ROOT types to python array type-codes.
    Below is a list of ROOT types:
            - C : a character string terminated by the 0 character
            - B : an 8 bit signed integer (Char_t)
            - b : an 8 bit unsigned integer (UChar_t)
            - S : a 16 bit signed integer (Short_t)
            - s : a 16 bit unsigned integer (UShort_t)
            - I : a 32 bit signed integer (Int_t)
            - i : a 32 bit unsigned integer (UInt_t)
            - F : a 32 bit floating point (Float_t)
            - D : a 64 bit floating point (Double_t)
            - L : a 64 bit signed integer (Long64_t)
            - l : a 64 bit unsigned integer (ULong64_t)
            - O : [the letter 'o', not a zero] a boolean (Bool_t)
    """
    conversion_dict = {'C':'c',
                       'B':'b',
                       'b':'B',
                       'S':'h',
                       's':'H',
                       'I':'i',
                       'i':'I',
                       'F':'f',
                       'D':'d',
                       'L':'l',
                       'l':'L',
                       'O':'b'}
    c_type = get_branch_ctype(branch)

    if c_type in conversion_dict:
        return conversion_dict[c_type]
    else:
        raise KeyError('rootutils.get_branch_pytype: Unknown conversion type')


def get_pyarray_initializer(pytype):
    """ returns the array initializer based on the python array type.
    Need this for SetBranchAddress call.
    """
    initializer_dict = {'c':'',
                        'b':[0],
                        'B':[0],
                        'h':[0],
                        'H':[0],
                        'i':[0],
                        'I':[0],
                        'f':[0.0],
                        'd':[0.0],
                        'l':[0],
                        'L':[0],
                        'b':[0]}

    return initializer_dict[pytype]


def sum_branch(chain, branch_name, branch_type=None, accessor=None):
    """Returns the sum of the entries of the branch 'branch_name' in TTree
    or TChain 'chain'

    branch_type: Optionally specify type of complex object for branch
    access

    accessor: Optionally specify accessor function for accessing values
    stored in complex type
    """
    if branch_type is None:
        branch_type = get_branch_pytype(chain.GetBranch(branch_name))
        branch_object = array(branch_type,
                               get_pyarray_initializer(branch_type))
        chain.SetBranchAddress(branch_name, branch_object)
    else:
        branch_object = branch_type()
        chain.SetBranchAddress(branch_name, rt.AddressOf(branch_object))

    total = 0
    for i in range(chain.GetEntries()):
        chain.GetEntry(i)

        if accessor is None:
            for a_value in branch_object:
                total += a_value
        else:
            total += accessor(branch_object)

    return total


def mean_branch(chain, branch_name):
    return sum_branch(chain, branch_name)/float(chain.GetEntries())


def atomZ_sel_str(pdg_var):
    """ Returns the ROOT sel string for cutting on Z using pdg codes.

    pdg_var: string representing the true-pdg variable in the TTree
    """
    return '{0}/10000 % 1000'.format(pdg_var)


def title(draw_func):
    """ decorator for turning on ROOT titles
    """
    @functools.wraps(draw_func)
    def wrap_func(*args, **kwargs):
        orig_top_margin = rt.gStyle.GetPadTopMargin()
        rt.gStyle.SetOptTitle(1)
        rt.gStyle.SetPadTopMargin(0.1)
        draw_func(*args, **kwargs)
        rt.gStyle.SetOptTitle(0)
        rt.gStyle.SetPadTopMargin(orig_top_margin)

    return wrap_func


def bigaxis(draw_func):
    """ decorator for increasing ROOT axis labels and title
    """
    @functools.wraps(draw_func)
    def wrap_func(*args, **kwargs):
        orig_title_size = rt.gStyle.GetTitleSize()
        orig_label_size = rt.gStyle.GetLabelSize()
        orig_x_offset = rt.gStyle.GetTitleOffset('x')
        orig_y_offset = rt.gStyle.GetTitleOffset('y')
        rt.gStyle.SetTitleSize(0.055, 'xy')
        rt.gStyle.SetLabelSize(0.05, 'xyz')
        rt.gStyle.SetTitleOffset(0.9, 'xy')
        draw_func(*args, **kwargs)
        rt.gStyle.SetTitleSize(orig_title_size, 'xy')
        rt.gStyle.SetLabelSize(orig_label_size, 'xyz')
        rt.gStyle.SetTitleOffset(orig_x_offset, 'x')
        rt.gStyle.SetTitleOffset(orig_y_offset, 'y')

    return wrap_func


def glob_to_chain(chain, paths):
    """ Default TChain::Add functionality doesn't support wildcards in
    directories, only base filename.  Workaround by adding individual
    files one by one.

    paths: wildcard or list of filepaths
    """
    if type(paths) is str:
        # glob wildcard
        paths = glob.glob(os.path.expanduser(paths))

    for path in paths:
        chain.AddFile(path)


def unix_time_axis_convert(taxis):
    """ Convert the given taxis labeling from unix timestamp to
    human readable dates.
    """
    taxis.SetTimeDisplay(1)
    taxis.SetTimeFormat("%Y/%m/%d%F1970-01-01 0:0:0")
    taxis.SetNdivisions(-5)


def binning(chain, branch, binstep=None):
    """ Sets the binning for histos to cover entire range of branch values in
    chain.

    Requires two arguments: chain and branch, which serve to find the max and
    min values of the branch.

    Returns an array of bin-low-X values that can be passed to TH1s
    """
    # protects against branch names that contain arrays
    branch = branch.split('[')[0]
    # automatically set the binstep if it's not passed as a parameter
    if binstep is None:
        binlow = 0.9*chain.GetMinimum(branch)
        binhigh = 1.1*chain.GetMaximum(branch)
        binstep = (binhigh-binlow)/50.
    else:
        binlow = chain.GetMinimum(branch)-binstep
        binhigh = chain.GetMaximum(branch)+binstep

    if binstep == 0:
        binstep = 1

    return array('f', utils.frange(binlow, binhigh+binstep, binstep))


def label_hist(hist, x_label='', y_label='', z_label='', tit=''):
    """ label root TH1 hist axes and titles.  The only required
    argument is 'hist'; all labels are empty by default.
    """
    hist.GetXaxis().SetTitle(x_label)
    hist.GetYaxis().SetTitle(y_label)
    hist.GetZaxis().SetTitle(z_label)
    hist.SetTitle(tit)


def inner_angle(v1, v2):
    """ v1, v2 are TVector3's.  Returns the inner angle between v1 and v2
    """
    return math.acos(v1.Dot(v2) / (v1.Mag() * v2.Mag()))


def split_varstring(varstring):
    """Splits a variable string containing branch names into a list with
    the single-colon as a delimiter.  Cannot simply use str.split()
    because functions like TMath::sin inject double colons into
    varstring.
    """
    # This regex checks that there are no colons before and after the
    # single colon
    return re.split('(?<!:):(?!:)', varstring)


def varstring_dimensions(varstring):
    """ Returns the dimensions of varstring (split on single colon)
    """
    return len(split_varstring(varstring))


def get_varx(varstring):
    """ Returns the x-axis variable in a ROOT varstring
    """
    split_vars = split_varstring(varstring)
    return (split_vars[1] if len(split_vars)==2
            else split_vars[0])


def get_vary(varstring):
    """ Returns the x-axis variable in a ROOT varstring
    For 1-D varstring return empty string
    """
    split_vars = split_varstring(varstring)
    # dim:vary dict
    vary_dict = {1:'',
                 2:split_vars[0],
                 3:split_vars[1],
                 4:split_vars[1]}

    return vary_dict[len(split_vars)]


def get_hist(name, hist_title, hist_binning):
    """tries to return the correct TH1 object based on the number of
    dimensions in hist_binning.  hist_binning is a tuple of arrays
    i.e. (binlowX, binlowY...)
    """
    dim = len(hist_binning)
    if dim == 1:
        return rt.TH1F(name, hist_title,
                       len(hist_binning[0])-1,
                       hist_binning[0])
    elif dim == 2:
        return rt.TH2F(name, hist_title,
                       len(hist_binning[0])-1,
                       hist_binning[0],
                       len(hist_binning[1])-1,
                       hist_binning[1])
    else:
        print ('ERROR! rootutils.get_hist: Not implmented for '
               'dimensions greater than 2!')


def dump_tree(tree, output, branch='*', sel=''):
    """Read and dump the tree entries to text file output.  Optionally
    specify a branch and selection
    """
    # use a copy so that we are not modifying state of tree and so can
    # still work with it interactively
    tree_copy = tree.Clone()

    tree_copy.GetPlayer().SetScanRedirect(True)
    tree_copy.GetPlayer().SetScanFileName(output)
    tree_copy.Scan(branch, sel)


def flatten_branch(orig_tree, branch_descriptor, branch_sel='*', sel=''):
    """This will flatten a branch using the TTree::Scan functionality,
    pipe the output to a txt file, and then read in the text file into
    a new flat tree.  This is useful if orig_tree's branches are
    object types or arrays and we want to expand them into own tree.

    branch_descriptor must contain a description of each branch in
    branch_sel.  It has syntax of the form:
    A/D:Table[2]/F:Ntracks/I:astring/C
    The 'Run' index and 'instance' (if it exists) will be automatically added
    """
    scan_out = os.path.join('~/.trash', 'flatten_branch_scan.txt')
    txt_out = os.path.join('~/.trash', 'flatten_branch_skim.txt')

    dump_tree(orig_tree, scan_out, branch_sel, sel)

    with open(scan_out) as scan_txt:
        out_txt = open(txt_out, 'w')

        # Remove ROOT Scan's header and asterisks which makes no sense to
        # ReadFile.
        for idx, line in enumerate(scan_txt):
            # Check for Instance column which should be in the second line
            if idx == 1:
                instance_descriptor = ('instance/I:'
                                       if line.split('*')[2].strip()=='Instance'
                                       else '')
            if idx > 2:
                out_txt.write(line.replace('*', ''))

        out_txt.close()

    branch_descriptor = '{0}{1}{2}'.format('orig_row/I:',
                                           instance_descriptor,
                                           branch_descriptor)
    new_tree = rt.TTree('new_tree', 'new tree')
    new_tree.ReadFile(txt_out, branch_descriptor,)

    os.remove(txt_out)
    os.remove(scan_out)
    return new_tree


def area_width_normalized_label(var, area_norm=True, width_norm=True):
    """Returns a label for histograms normalized by any combination of
    area/width. For both area and width, the label should be 1/N*dN/d(var).
    """
    if area_norm and width_norm:
        return '#frac{dN}{Nd('+var+')}'
    elif area_norm:
        return '#frac{#Delta('+var+')dN}{Nd('+var+')}'
    elif width_norm:
        return '#frac{dN}{d('+var+')}'
    else:
        return '#frac{#Delta('+var+')dN}{d('+var+')}'


def escape_for_filename(rootstr):
    """Given a selection rootstr string, returns a properly formatted string
    for filenaming.
    """
    if rootstr.strip()=='':
        return 'none'

    return (rootstr.
            lower().
            replace(' ', '').
            replace('>', 'gt').
            replace('<', 'lt').
            replace('==', 'eq').
            replace('&&', 'a').
            replace('||', 'o').
            replace('!', 'n').
            replace('[', '').
            replace(']', '').
            replace('(', 'p').
            replace(')', 'p').
            replace('_', '-'))


def get_line(const, canv, color=rt.kGray, width=1, style=9, vertical=False):
    """Gets a (default horizontal) TLine at const that stretches across
    the entire range of canv
    """
    # must call update to ensure correct canvas ranges are set by TH1::Draw
    canv.Update()

    if vertical:
        return rt.TLine(const, canv.GetUymin(), const, canv.GetUymax())
    else:
        return rt.TLine(canv.GetUxmin(), const, canv.GetUxmax(), const)


def zeroline(canv):
    """ returns a dashed, gray horizontal line at zero
    """
    zline = get_line(0, canv)
    zline.SetLineColor(rt.kGray)
    zline.SetLineStyle(9)
    return zline


def bin_log(nbins, bin_min, bin_max, non0start=0.01, force_min=False):
    """Sets normal_bins to log binning. This makes the bin widths appear
    the same length on a log-scale plot. Ported from
    http://www-pnp.physics.ox.ac.uk/~luxi/transport/Codes/2015_02_0607/T2K/NeutrinoTools.h

    non0start: forces the log binning to start at a set value greater
    than 0. If set to 0, force the use of bin_min as the start of the
    binning range by adding it on in front. Beware that this will make
    it impossible to draw equal width log-scale bins.
    """
    bin_min_orig = bin_min
    # avoid trying to draw 0 on log-scale
    if bin_min < sys.float_info.epsilon:
        bin_min = non0start

    factor = pow(bin_max/bin_min, 1./nbins)

    log_bins = [bin_min*pow(factor, i) for i in range(nbins+1)]
    return [bin_min_orig]+log_bins if force_min else log_bins


def tgt_color_code(number):
    """ Returns a root color for tgt with atomic number 'number'
    """
    tgt_color_dict = {6:rt.kRed,
                      8:rt.kBlue,
                      82:rt.kOrange}

    return tgt_color_dict[number] if number in tgt_color_dict else number


def get_flat1d_bin(bin_nd, hist, inc_overflow):
    """Converts global bin number bin_nd to the corresponding flattened
    1d global bin number. Can pass option to include nd overflow bins
    as individual bins in 1d
    """
    nbinsx = hist.GetNbinsX()
    nbinsy = hist.GetNbinsY()
    nbinsz = hist.GetNbinsZ()

    if inc_overflow:
        # Must +1 here to ensure that bin_nd_0 -> bin_1d_1 since it's
        # converted to it's own bin
        return bin_nd+1
    else:
        # Get the corresponding x, y, z to bin_nd
        binx = rt.Long(0)
        biny = rt.Long(0)
        binz = rt.Long(0)
        hist.GetBinXYZ(bin_nd, binx, biny, binz)

        if ((binx == 0 and nbinsx > 1) or
            (biny == 0 and nbinsy > 1) or
            (binz == 0 and nbinsz > 1)):
            # underflow
            return 0
        elif (binx > nbinsx > 1 or
              biny > nbinsy > 1 or
              binz > nbinsz > 1):
            # overflow
            return nbinsx*nbinsy*nbinsz + 1
        else:
            if nbinsy == 1:
                return binx
            elif nbinsz == 1:
                return binx + nbinsx*(biny-1)
            else:
                return binx + nbinsx*(biny-1) + nbinsx*nbinsy*(binz-1)


def get_flat1d_size(hist, inc_overflow):
    """ Returns the appropriate size of the flattend 1d histo
    """
    nbinsx = hist.GetNbinsX()
    nbinsy = hist.GetNbinsY()
    nbinsz = hist.GetNbinsZ()

    return nbinsx*nbinsy*nbinsz if not inc_overflow else hist.GetSize()


def is_underflow(bin_nd, hist):
    """Retuns whether global bin number bin_nd is an underflow bin. Works
    for any number of dimensions
    """
    flat1d_bin = get_flat1d_bin(bin_nd, hist, False)
    return flat1d_bin == 0


def is_overflow(bin_nd, hist):
    """Retuns whether global bin number bin_nd is an overflow bin. Works
    for any number of dimensions
    """
    flat1d_bin = get_flat1d_bin(bin_nd, hist, False)
    return flat1d_bin == get_flat1d_size(hist, False)+1


def flatten_hist(hist, inc_overflow=False):
    """Flattens a nD histogram into a 1D histo

    if inc_overflow is True, each overflow bin will be treated as a
    single 1D bin, otherwise it will be pushed to the overflow and
    underflow 1D bins
    """
    nbins1d = get_flat1d_size(hist, inc_overflow)
    binlow = 0 if inc_overflow else 1
    binhigh = nbins1d if inc_overflow else nbins1d+1

    hflat = rt.TH1D('hflat',
                    hist.GetTitle(),
                    nbins1d, binlow, binhigh)

    # loop over all bins in the nd hist
    for bin_nd in range(hist.GetSize()):
        bin_1d = get_flat1d_bin(bin_nd, hist, inc_overflow)

        ## testing
        # binx = rt.Long(0)
        # biny = rt.Long(0)
        # binz = rt.Long(0)
        # hist.GetBinXYZ(bin_nd, binx, biny, binz)
        # print binx, biny, bin_1d
        ## testing

        hflat.SetBinContent(bin_1d,
                            hflat.GetBinContent(bin_1d)+hist.GetBinContent(bin_nd))
        hflat.SetBinError(bin_1d,
                          math.sqrt(
                              hflat.GetBinError(bin_1d)**2+hist.GetBinError(bin_nd)**2))


    return hflat


def projection_normalize(orig_hist, xaxis=True):
    """Normalize all bins across the projection relative to the
    projection's max
    """
    hist = orig_hist.Clone()
    nbinsx = hist.GetNbinsX()
    nbinsy = hist.GetNbinsY()

    walkaxis = nbinsy if xaxis else nbinsx
    projaxis = nbinsx if xaxis else nbinsy

    # walk over bins in non-projected axis
    for b in range(walkaxis):
        projected = (hist.ProjectionX('_px', b, b+1) if xaxis else
                     hist.ProjectionY('_py', b, b+1))
        maxprojected = max(abs(projected.GetMaximum()),
                           abs(projected.GetMinimum()))

        if maxprojected == 0:
            continue

        # renormalize by stepping over bins in axis we project on
        for i in range(projaxis):
            xbin, ybin = (i, b) if xaxis else (b, i)
            hist.SetBinContent(xbin, ybin,
                               hist.GetBinContent(xbin, ybin)/maxprojected)
            hist.SetBinError(xbin, ybin,
                             hist.GetBinError(xbin, ybin)/maxprojected)

    return hist


def hist_to_formula(hist, var, sel='1'):
    """Converts a reweighting histogram into a TTree formula so that the
    formula can be used as a selection (e.g. in TTree::Draw) instead
    of going through event-by-event and reweighting using the
    histogram.

    parameters:
        hist- histogram to be used for reweighting
        var- variable along the histogram's x-axis (currently only 1d supported)
        sel- additional selection criteria
    returns: formula
    """
    # increase the maximum number of operations a formula can have (default 1000)
    rt.TFormula.SetMaxima(20000)

    bin_formulas = []
    # include under/overflow
    xaxis = hist.GetXaxis()
    for b in range(hist.GetNbinsX()+2):
        low = xaxis.GetBinLowEdge(b)
        up = xaxis.GetBinUpEdge(b)
        w = hist.GetBinContent(b)
        if b==0:
            bin_formulas.append('{0}*({1} < {2})'.format(w, var, up))
        elif b==hist.GetNbinsX()+1:
            bin_formulas.append('{0}*({1} >= {2})'.format(w, var, low))
        else:
            bin_formulas.append('{0}*({1} >= {2} && {1} < {3})'.format(w, var, low, up))

    return '(({0})*({1})+1*!({0}))'.format(sel, '+'.join(bin_formulas))


def contents(hist):
    """ returns the histograms contents in a python list
    """
    return [hist.GetBinContent(i) for i in range(hist.GetSize())]


def errors(hist):
    """ returns the histograms contents in a python list
    """
    return [hist.GetBinError(i) for i in range(hist.GetSize())]


def edges(hist):
    """ returns the histograms contents in a python list
    """
    return [hist.GetBinLowEdge(i) for i in range(hist.GetSize()+1)]


def tabulate(hist):
    """returns a list of tuples that represents the histogram's data in
    the format [(bin, content, error),...]
    """
    return list(zip(list(range(hist.GetSize())), contents(hist), errors(hist)))


def chi2test(hobs, hexp, opts=''):
    """Alternative to root's chi2 test. This takes the normalization
    differences into account and acts as an absolute chi2 test. Errors
    are added in quadrature. zero bins are skipped.
    """
    chi2 = 0
    ndof = 0
    for b in range(hobs.GetSize()):
        # skip underflow/overflow bins by default
        if is_underflow(b, hobs)  and 'UF' not in opts:
            continue
        if is_overflow(b, hobs) and 'OF' not in opts:
            continue
        obs_c = hobs.GetBinContent(b)
        exp_c = hexp.GetBinContent(b)
        obs_e = hobs.GetBinError(b)
        exp_e = hexp.GetBinError(b)

        if obs_c == 0 and exp_c == 0:
            continue
        if obs_e == 0 and exp_e == 0:
            continue
        chi2 += (obs_c - exp_c)**2/(obs_e**2+exp_e**2)
        ndof += 1

    pval = rt.TMath.Prob(chi2, ndof)

    if 'p' in opts.lower():
        print('INFO: Chi2 = {0}, Prob = {1}, NDF = {2}'.format(chi2,
                                                               pval,
                                                               ndof))

    if 'chi2/ndf' in opts.lower():
        return chi2/ndof
    elif 'chi2' in opts.lower():
        return chi2
    else:
        return pval


def reduce_cvm(cvm):
    """ removes rows and columns along that correspond to 0s on the diagonal
    """
    nrows = cvm.GetNrows()
    assert nrows == cvm.GetNcols()
    diag = rt.TMatrixDDiag(cvm)
    nnonzeros = 0
    inonzeros = []
    for i in range(diag.GetNdiags()):
        if diag[i] != 0:
            nnonzeros += 1
            inonzeros.append(i)

    to_save = []
    for i in inonzeros:
        for j in inonzeros:
            to_save.append(cvm[i][j])

    rcvm = rt.TMatrixD(nnonzeros, nnonzeros)
    overall = 0
    for i in range(nnonzeros):
        for j in range(nnonzeros):
            rcvm[i][j] = to_save[overall]
            overall += 1

    return rcvm
