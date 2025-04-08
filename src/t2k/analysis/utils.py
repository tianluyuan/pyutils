import os
import re
import ROOT as rt
from tianlu import directories, rootutils
from tianlu import utils as tutils
from tianlu.analysis.varstrings import *
from tianlu.analysis.binning import *
from distutils.spawn import find_executable
from collections import namedtuple
import ordereddict


# create a namedtuple to pair up reco and true cuts
SELECTIONPAIR = namedtuple('SELECTIONPAIR', 'reco true')


def nd280analysis_prod():
    """ Get the production version used in compiling nd280AnalysisTools
    """
    nd280_input_file = os.path.join(os.getenv('ND280ANALYSISTOOLSROOT'),
                                    'AnalysisTools/input-file.list')

    with open(nd280_input_file) as f:
        oaanal_file = os.path.basename(f.readline().strip('/n'))

    return int(re.findall('prod([0-9]*)', oaanal_file)[0])


def flux_ratio_file():
    """ Returns the appropriate file for flux reweighting based on
    the oaAnalysis version used to compile nd280AnalysisTools.
    Flux ratios are different for prod6 and we want to enforce that
    analysis scripts such as runAnalyses.py -a GetEvents use the
    correct ratios depending on the production being analyzed.

    This should probably be made more robust.
    """
    prod_version = nd280analysis_prod()

    if prod_version == 5:
        # prod5 nominal flux is 11a
        return os.path.join(directories.NFS_BASE,
                            '/beam/flux/13a/tuned13av1.0/run1-4/nd5_tuned13av1.0_11anom_run1-4_fine.root')
    elif prod_version == 6:
        # prod6 nominal flux is 13a
        return os.path.join(directories.NFS_BASE,
                            '/beam/flux/13a/tuned13av1.0/run1-4/nd5_tuned13av1.0_13anom_run1-4_fine.root')
    else:
        print 'Undefined flux ratios file for production', prod_version
        print 'Using 13a/11b ratios'
        return os.path.join(os.getenv('BASEANALYSISROOT'),
                            'data/nd5_11btune_13anom_250ka_ratios.root')


def env_highland_version():
    """Get the highland version currently loaded in the environment by
    accessing HIGHLANDTOOLSROOT.  This is useful when we're running a
    program that needs highland files as input but does not have a
    path that tells us which highland version being used.
    """
    return tutils.split_path(os.getenv('HIGHLANDTOOLSROOT'))[-3]


def env_is_highland1():
    return env_highland_version() == 'highland'


def env_is_highland2():
    return env_highland_version() == 'highland2'


def exe_highland_version(exe):
    """ Checks if the executable program is highland2 or highland.  If
    neither throw an exception.
    """
    path = find_executable(exe)
    if 'highland2' in path:
        return 2
    elif 'highland' in path:
        return 1
    else:
        raise RuntimeError(exe+' is not a highland program!')


def xsthrow_format(formula):
    """formats the string to follow the xstool_throw convention for toy
    vars
    """
    return (formula.
            replace('accum_level[0]', 'accum_level[xstool_throw]').
            replace('selmu_mom[0]', 'selmu_mom[xstool_throw]').
            replace('selmu_theta[0]', 'selmu_theta[xstool_throw]'))


def convert_var_name(var):
    """Returns a nice variable name as oppposed to the ROOT tree variable
    itself for labeling
    """
    highland_var_dict = {TRUE_NU_E:'E_{#nu} [GeV]',
                         TRUE_MU_MOM:'True-#mu p_{#mu} [GeV]',
                         TRUE_MU_E:'True-#mu E_{#mu} [GeV]',
                         TRUE_MU_POSZ:'True-#mu z_{#mu} [mm]',
                         TRUE_MU_CTH:'True-#mu cos#theta_{#mu}',
                         TRUE_MU_TH:'True-#mu #theta_{#mu}',
                         TRUE_Q2:'Q^{2} [GeV^{2}]',
                         TRUE_X:'x_{Bj}',
                         TRUE_P_MOMT:'True-p p_{T} [GeV]',
                         TRUE_MU_MOMT:'True-#mu p_{T} [GeV]',
                         TRUE_DPHI:'#delta#phi_{T}',
                         SEL_DPHI:'Reco-#delta#phi_{T}',
                         SEL_MU_MOMT:'Reco-#mu p_{T} [GeV]',
                         SEL_MU_MOM:'Reco-#mu p_{#mu} [GeV]',
                         SEL_MU_CTH:'Reco-#mu cos#theta_{#mu}',
                         SEL_MU_TH:'Reco-#mu #theta_{#mu}',
                         SEL_MU_POSZ:'Reco-#mu z_{#mu} [mm]',
                         DELTA_DPHI:'Reco-#delta#phi_{T}-True-#delta#phi_{T}'}

    # apply the xsthrow_format function across all keys in dict and
    # build a new var_dict containing the full key value pair
    var_dict = {}
    for key, value in highland_var_dict.iteritems():
        var_dict[key] = value
        var_dict[xsthrow_format(key)] = value

    if var_dict.has_key(var):
        return var_dict[var]
    else:
        return var


def xsec_label_from_varstring(varstring):
    """sets the xsec label (e.g. dsigma/dp) based on the varstring.
    Defined for all variables besides TRUE_NU_E, which is
    flux-normalized and hence just sigma(E).  The dictionary contains
    tuple-values with the first entry being the variables abbreviation
    and the second being the units.
    """
    xsec_label_dict = {TRUE_MU_MOM:('p','/GeV'),
                       TRUE_MU_MOMT:('p_{T}','/GeV'),
                       TRUE_MU_E:('E_{#mu}','/GeV'),
                       TRUE_MU_POSZ:('z', '/mm'),
                       TRUE_MU_CTH:('cos#theta', ''),
                       TRUE_MU_TH:('#theta', ''),
                       TRUE_DPHI:('#delta#phi_{T}', ''),
                       TRUE_Q2:('Q^{2}', '/GeV^{2}'),
                       TRUE_X:('x_{Bj}', '')}

    varlist = rootutils.split_varstring(varstring)
    dim = len(varlist)
    if dim == 1:
        varabbrv = xsec_label_dict[varlist[0]][0]
        varunits = xsec_label_dict[varlist[0]][1]
        label = ('#frac{d#sigma}{d'+varabbrv+'} '
                '(10^{-38} cm^{2}'+varunits+'/neutron)')
    elif dim == 2:
        varabbrv = (xsec_label_dict[varlist[0]][0],
                    xsec_label_dict[varlist[1]][0])
        varunits = (xsec_label_dict[varlist[0]][1],
                    xsec_label_dict[varlist[1]][1])
        label = ('#frac{d^{2}#sigma}{d'
                 ''+varabbrv[0]+'d'
                 ''+varabbrv[1]+'} (10^{-38} cm^{2}'
                 ''+varunits[0]+''
                 ''+varunits[1]+'/neutron)')

    return label


def set_axes_titles(varstring, histo, title=''):
    """Labels axes based on the varstring.  Tries to be smart about the
    number of dimensions and sets axes appropriately. Converts
    varstring to a nice readable format.  The last arg 'title' is
    optional and defaults to empty string.
    """
    # use regex to match varstrings for two dimensions'
    dim = rootutils.varstring_dimensions(varstring)
    split_vars = rootutils.split_varstring(varstring)

    if dim == 1:
        histo.GetXaxis().SetTitle(
            convert_var_name(split_vars[0]))
        histo.GetYaxis().SetTitle(title)
    elif dim == 2:
        histo.GetXaxis().SetTitle(
            convert_var_name(split_vars[1]))
        histo.GetYaxis().SetTitle(
            convert_var_name(split_vars[0]))
        histo.GetZaxis().SetTitle(title)
    elif dim == 3:
        histo.GetXaxis().SetTitle(
            convert_var_name(split_vars[0]))
        histo.GetYaxis().SetTitle(
            convert_var_name(split_vars[1]))
        histo.GetZaxis().SetTitle(
            convert_var_name(split_vars[2]))

    histo.SetTitle(title)


def prettify_highland_filepath(filepath):
    """ Make highland output filename prettier
    RunP0DNumuCCQEAnalysis_prod6B_mcp_neut_RUN1_merged.root -> 'NEUT RUN1'
    """
    tokenized = os.path.splitext(os.path.basename(filepath))[0].split('_')

    return ' '.join(tokenized[3:5]).upper()


def draw_binning(binning=MOM_CTH_PI2_2D,
                 labels=(r'$p_\mu$ (GeV/c)', r'$\cos \theta_\mu$'),
                 gbins=False):
    """ Draws a 2D-binning in a nice viewable format.

    binning: A 2D-tuple of array-binnings in form (x-edges, y-edges)
    labels: A 2D-tuple of axis labels
    gbins: Include global bin numbering
    """
    # import here so that this utils file can be imported even if
    # matplotlib isn't installed, as matplotlib isn't always available
    # on all machines
    import matplotlib.pyplot as plt

    if len(binning) != 2:
        raise RuntimeError('Not implemented for dimensions != 2')

    plt.rc('text', usetex=True)
    plt.xticks(binning[0], size=20)
    plt.yticks(binning[1], size=20)
    plt.xlabel(labels[0], size=20)
    plt.ylabel(labels[1], size=20)
    plt.grid(linestyle='-', linewidth=2)
    if gbins:
        midx = tutils.midpoints(binning[0])
        midy = tutils.midpoints(binning[1])
        for idy, y in enumerate(midy):
            for idx, x in enumerate(midx):
                plt.text(0.98*x, 0.99*y, idy*len(midx) + idx + 1,
                         color='r', weight='heavy', size=14)
    plt.show()


def unstack_sources(stacked):
    """ Unstacks the stacked sources list, where stacked sources are given
    as a ordereddict of tuples in the form set in highland_configurations
    """
    categories = stacked.keys()

    unstacked = ordereddict.OrderedDict()

    for idx, key in enumerate(categories[:-1]):
        # take the difference between two stacked groups to be the
        # unstacked list of sources for that particular category
        nkey = categories[idx+1]
        unstacked[key] = (list(set(stacked[key][0])-set(stacked[nkey][0])),
                          stacked[key][1])

    # the last category in the ordered dict is its own unstacked group
    lkey = categories[-1]
    unstacked[lkey] = stacked[lkey]

    return unstacked


def get_nominal(param, banff_file, antinumode=False, postfit=False):
    """ returns the nominal value of a parameter as set in xsReadBanffFile
    """
    rbf = rt.xsReadBanffFile(banff_file, antinumode, postfit)
    dials = list(rbf.GetDialNames())
    vals = list(rbf.GetNominalDialValues())
    return vals[dials.index(param)]


def get_variance(param, banff_file, antinumode=False, postfit=False):
    """ returns the nominal value of a parameter as set in xsReadBanffFile
    """
    rbf = rt.xsReadBanffFile(banff_file, antinumode, postfit)
    dials = list(rbf.GetDialNames())
    cvm = rbf.GetCovarianceMatrix()
    idx = dials.index(param)
    return cvm[idx][idx]
