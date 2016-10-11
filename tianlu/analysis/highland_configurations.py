"""
Options file that can be imported and used for unfolding.  Good to separate
data from behavior.
"""
import os
import ROOT as rt
import ordereddict
from tianlu import directories
from tianlu.analysis.utils import SELECTIONPAIR
from tianlu.analysis.P0DCutUtils import getLayerCutStr
from collections import defaultdict


__author__ = 'Tianlu Yuan'
__email__ = 'tianlu.yuan [at] colorado.edu'


#################### Globals ######################
# Set debug to false when analysis is finalized and data can be used
DEBUG = False
# Use T2K style for DrawingTools
T2KSTYLE = False
# production version to use as inputs
PROD_MCP = 'prod6B'
PROD_RDP = 'prod6E'
# highland analysis specifier. need to change SIGNAL and CUTS if this changes
ANALYSIS = 'RunP0DNumuCCQEAnalysis'
# Set to True to use merged MC files, otherwise will use a single
# output instead of merged. If true, a cache file is saved to speed up
# later processing if xsToolsUnfold.stack_sources is called
FULLSTATS = True
# flag for whether to set the fake-data to data POT
FAKEPOT = True
#################### Inputs ######################
FILE_SPEC = 'merged' if FULLSTATS else '1'

# MCP NEUT
NEUT_TEMPLATE = os.path.join(directories.HIGHLAND2_STORE, ANALYSIS, PROD_MCP,
                             'mcp', 'neut', '{0}',
                             '{0}_{1}_{2}'.format(ANALYSIS,
                                                  PROD_MCP,
                                                  'mcp_neut_{0}_'+FILE_SPEC+'.root'))
MCP_NEUT_RUN1 = NEUT_TEMPLATE.format('RUN1')
MCP_NEUT_RUN2A = NEUT_TEMPLATE.format('RUN2air')
MCP_NEUT_RUN2W = NEUT_TEMPLATE.format('RUN2water')
MCP_NEUT_RUN3C = NEUT_TEMPLATE.format('RUN3')
MCP_NEUT_RUN4W = NEUT_TEMPLATE.format('RUN4water')
MCP_NEUT_RUN4A = NEUT_TEMPLATE.format('RUN4air')

# MCP GENIE
GENIE_TEMPLATE = os.path.join(directories.HIGHLAND2_STORE, ANALYSIS, PROD_MCP,
                              'mcp', 'genie', '{0}',
                              '{0}_{1}_{2}'.format(ANALYSIS,
                                                   PROD_MCP,
                                                   'mcp_genie_{0}_'+FILE_SPEC+'.root'))
MCP_GENIE_RUN1 = GENIE_TEMPLATE.format('RUN1')
MCP_GENIE_RUN2A = GENIE_TEMPLATE.format('RUN2air')
MCP_GENIE_RUN2W = GENIE_TEMPLATE.format('RUN2water')
MCP_GENIE_RUN3C = GENIE_TEMPLATE.format('RUN3')
MCP_GENIE_RUN4W = GENIE_TEMPLATE.format('RUN4water')
MCP_GENIE_RUN4A = GENIE_TEMPLATE.format('RUN4air')

# RDP
BEAM_TEMPLATE = os.path.join(directories.HIGHLAND2_STORE, ANALYSIS, PROD_RDP,
                             'rdp', 'beam', '{0}',
                             '{0}_{1}_{2}'.format(ANALYSIS,
                                                  PROD_RDP,
                                                  'rdp_beam_{0}_'+FILE_SPEC+'.root'))
RDP_RUN1 = BEAM_TEMPLATE.format('RUN1')
RDP_RUN2A = BEAM_TEMPLATE.format('RUN2air')
RDP_RUN2W = BEAM_TEMPLATE.format('RUN2water')
RDP_RUN3C = BEAM_TEMPLATE.format('RUN3c')
RDP_RUN4W = BEAM_TEMPLATE.format('RUN4water')
RDP_RUN4A = BEAM_TEMPLATE.format('RUN4air')

# Weights
# t2krw
WEIGHTS_TEMPLATE = os.path.join(ANALYSIS, PROD_MCP, 'mcp', 'neut',
                                '{0}', FILE_SPEC, 'weights_unified.root')
WEIGHTS_NEUT_RUN1 = WEIGHTS_TEMPLATE.format('RUN1')
WEIGHTS_NEUT_RUN2W = WEIGHTS_TEMPLATE.format('RUN2water')
WEIGHTS_NEUT_RUN2A = WEIGHTS_TEMPLATE.format('RUN2air')
WEIGHTS_NEUT_RUN3C = WEIGHTS_TEMPLATE.format('RUN3')
WEIGHTS_NEUT_RUN4W = WEIGHTS_TEMPLATE.format('RUN4water')
WEIGHTS_NEUT_RUN4A = WEIGHTS_TEMPLATE.format('RUN4air')
# mass
MASS_TEMPLATE = os.path.join(ANALYSIS, PROD_MCP, 'mcp', 'neut',
                             '{0}', FILE_SPEC, 'weights_mass.root')
MASS_NEUT_RUN1 = MASS_TEMPLATE.format('RUN1')
MASS_NEUT_RUN2W = MASS_TEMPLATE.format('RUN2water')
MASS_NEUT_RUN2A = MASS_TEMPLATE.format('RUN2air')
MASS_NEUT_RUN3C = MASS_TEMPLATE.format('RUN3')
MASS_NEUT_RUN4W = MASS_TEMPLATE.format('RUN4water')
MASS_NEUT_RUN4A = MASS_TEMPLATE.format('RUN4air')

# xsTool generated rwfiles, needed to get NRooTrackerVtxs
RWFILE_DEFAULT_TEMPLATE = os.path.join(ANALYSIS, PROD_MCP, 'mcp', 'neut',
                                       '{0}', FILE_SPEC, 'rwfile_default.root')
RWFILE_DEFAULT_NEUT_RUN1 = RWFILE_DEFAULT_TEMPLATE.format('RUN1')
RWFILE_DEFAULT_NEUT_RUN2W = RWFILE_DEFAULT_TEMPLATE.format('RUN2water')
RWFILE_DEFAULT_NEUT_RUN2A = RWFILE_DEFAULT_TEMPLATE.format('RUN2air')
RWFILE_DEFAULT_NEUT_RUN3C = RWFILE_DEFAULT_TEMPLATE.format('RUN3')
RWFILE_DEFAULT_NEUT_RUN4W = RWFILE_DEFAULT_TEMPLATE.format('RUN4water')
RWFILE_DEFAULT_NEUT_RUN4A = RWFILE_DEFAULT_TEMPLATE.format('RUN4air')

RWFILE_TRUTH_TEMPLATE = os.path.join(ANALYSIS, PROD_MCP, 'mcp', 'neut',
                                     '{0}', FILE_SPEC, 'rwfile_truth.root')
RWFILE_TRUTH_NEUT_RUN1 = RWFILE_TRUTH_TEMPLATE.format('RUN1')
RWFILE_TRUTH_NEUT_RUN2W = RWFILE_TRUTH_TEMPLATE.format('RUN2water')
RWFILE_TRUTH_NEUT_RUN2A = RWFILE_TRUTH_TEMPLATE.format('RUN2air')
RWFILE_TRUTH_NEUT_RUN3C = RWFILE_TRUTH_TEMPLATE.format('RUN3')
RWFILE_TRUTH_NEUT_RUN4W = RWFILE_TRUTH_TEMPLATE.format('RUN4water')
RWFILE_TRUTH_NEUT_RUN4A = RWFILE_TRUTH_TEMPLATE.format('RUN4air')
##################### End #######################

##################### Flux ######################
FLUX_RUN1 = '{0}/beam/flux/13a/tuned13av1.1/run1/nd5_tuned13av1.1_13anom_run1_fine.root'.format(directories.NFS_BASE)
FLUX_RUN2 = '{0}/beam/flux/13a/tuned13av1.1/run2/nd5_tuned13av1.1_13anom_run2_fine.root'.format(directories.NFS_BASE)
FLUX_RUN3C = '{0}/beam/flux/13a/tuned13av1.1/run3c/nd5_tuned13av1.1_13anom_run3c_fine.root'.format(directories.NFS_BASE)
FLUX_RUN4 = '{0}/beam/flux/13a/tuned13av1.1/run4/nd5_tuned13av1.1_13anom_run4_fine.root'.format(directories.NFS_BASE)
FLUX_HIST_NUMU = 'enu_nd5_tuned13a_numu'
# for xsReweightSystematics.  This specifies which neutrino is used for the
# flux weights
NEUTRINO = 'numu'

################ fake-data POT ##################
# for fake-data studies
# taken from the actual highland2 merged data outputs and retrieved
# using xsInputdataHighland2::GetPOTdata
if FAKEPOT:
     FAKEPOT_RUN1 = 1.67e19
     FAKEPOT_RUN2W = 4.29e19
     FAKEPOT_RUN2A = 3.55e19
     FAKEPOT_RUN3C = 1.35e20
     FAKEPOT_RUN4W = 1.63e20
     FAKEPOT_RUN4A = 1.76e20
else:
    FAKEPOT_RUN1 = -1
    FAKEPOT_RUN2W = -1
    FAKEPOT_RUN2A = -1
    FAKEPOT_RUN3C = -1
    FAKEPOT_RUN4W = -1
    FAKEPOT_RUN4A = -1
##################### End #######################

################## Selection ####################
# (T.Y. 5/6/14) && particle ==13'
SELECT = SELECTIONPAIR('Sum$(accum_level[0][0]) > 3', 'topology==0')
INCLUSIVE = SELECTIONPAIR('Sum$(accum_level[0][0]) > 2', 'topology<3')
SELECT_2TRACK0MICHEL = SELECTIONPAIR('Sum$(accum_level[0][1]) > 3 && nmichel==0 && np0d==2', 'topology==0')
##################### End #######################

################### Sidebands ####################
# list all sidebands to be added here
SIDEBANDS = [SELECTIONPAIR('Sum$(accum_level[0][1]) > 3 && np0d == 2 && nmichel>0', 'topology == 1'),
             SELECTIONPAIR('Sum$(accum_level[0][1]) > 3 && np0d > 2', 'topology == 2')]
# disable for specific sources
SIDEBAND_DISABLED = ['T2KReWeight_flux', 'Mass_all']
##################### End #######################

################# Backgrounds ###################
# Only the background components.
SELECT_BACKGROUND = SELECTIONPAIR('Sum$(accum_level[0][0]) > 3 && topology!=0',
                                  'topology!=0')
SELECT_CC1PI = SELECTIONPAIR('Sum$(accum_level[0][0]) > 3 && topology==1',
                             'topology==1')
SELECT_CCOTHER = SELECTIONPAIR('Sum$(accum_level[0][0]) > 3 && topology==2',
                               'topology==2')
SIDEBAND1_CC1PI = SELECTIONPAIR('&&'.join(SIDEBANDS[0]),
                                SIDEBANDS[0].true)
SIDEBAND2_CCOTHER = SELECTIONPAIR('&&'.join(SIDEBANDS[1]),
                                  SIDEBANDS[1].true)
##################### End #######################

################ Inefficiency ###################
MISSED = SELECTIONPAIR('Sum$(accum_level[0][0]) <=3 && topology==0', 'topology==0')
##################### End #######################

########### Categories Sideband Fit #############
# for category fitting
CATEGS = ['topology==1', 'topology==2']
##################### End #######################

##################### Errors ####################
COLORS = defaultdict(lambda:rt.kOrange,
                     Mass_all=rt.kOrange+2,
                     statistics_data=rt.kCyan,
                     statistics_mc=rt.kCyan+2,
                     T2KReWeight_fsi=rt.kRed,
                     T2KReWeight_xs=rt.kYellow,
                     T2KReWeight_flux=rt.kGreen)
# List here the groups of error sources that we want to include.
# Sort in order of decreasing number of sources for drawing on top of previous
# Empty list corresponds to all implemented errors.
MASTER = ordereddict.OrderedDict()
MASTER['detector'] = (['statistics_data',
                       'statistics_mc',
                       'T2KReWeight_fsi',
                       'T2KReWeight_xs',
                       'T2KReWeight_flux',
                       'Mass_all',
                       'highland_all_syst'], COLORS['detector'])
MASTER['mass'] = (['statistics_data',
                   'statistics_mc',
                   'T2KReWeight_fsi',
                   'T2KReWeight_xs',
                   'T2KReWeight_flux',
                   'Mass_all'], COLORS['Mass_all'])
MASTER['flux'] = (['statistics_data',
                   'statistics_mc',
                   'T2KReWeight_fsi',
                   'T2KReWeight_xs',
                   'T2KReWeight_flux'], COLORS['T2KReWeight_flux'])
MASTER['cross-section'] = (['statistics_data',
                            'statistics_mc',
                            'T2KReWeight_fsi',
                            'T2KReWeight_xs'], COLORS['T2KReWeight_xs'])
MASTER['fsi'] = (['statistics_data',
                  'statistics_mc',
                  'T2KReWeight_fsi'], COLORS['T2KReWeight_fsi'])
MASTER['statistical (mc)'] = (['statistics_data',
                               'statistics_mc'], COLORS['statistics_mc'])
MASTER['statistical (data)'] = (['statistics_data'], COLORS['statistics_data'])

STATSYST = ordereddict.OrderedDict()
STATSYST['systematics'] = (['statistics_data',
                            'statistics_mc',
                            'T2KReWeight_fsi',
                            'T2KReWeight_xs',
                            'T2KReWeight_flux',
                            'Mass_all',
                            'highland_all_syst'], COLORS['T2KReWeight_fsi'])
STATSYST['statistics'] = (['statistics_data',
                           'statistics_mc'], COLORS['statistics_data'])

STATONLY = ordereddict.OrderedDict()
STATONLY['statistical (mc)'] = (['statistics_data',
                                 'statistics_mc'], COLORS['statistics_mc'])
STATONLY['statistical (data)'] = (['statistics_data'], COLORS['statistics_data'])

STATDATA = ordereddict.OrderedDict()
STATDATA['statistical (data)'] = (['statistics_data'], COLORS['statistics_data'])

STATMC = ordereddict.OrderedDict()
STATMC['statistical (mc)'] = (['statistics_mc'], COLORS['statistics_mc'])

XSFSIMC = ordereddict.OrderedDict()
XSFSIMC['fsi'] = (['statistics_mc',
                   'T2KReWeight_xs',
                   'T2KReWeight_fsi'], COLORS['T2KReWeight_fsi'])
XSFSIMC['cross-section'] = (['statistics_mc',
                             'T2KReWeight_xs'], COLORS['T2KReWeight_xs'])
XSFSIMC['statistical (mc)'] = (['statistics_mc'], COLORS['statistics_mc'])

MODELSTAT = ordereddict.OrderedDict()
MODELSTAT['cross-section'] = (['statistics_data',
                               'statistics_mc',
                               'T2KReWeight_xs'], COLORS['T2KReWeight_xs'])
MODELSTAT['statistical (mc)'] = (['statistics_data',
                                 'statistics_mc'], COLORS['statistics_mc'])
MODELSTAT['statistical (data)'] = (['statistics_data'], COLORS['statistics_data'])

FLUXSTAT = ordereddict.OrderedDict()
FLUXSTAT['flux'] = (['statistics_data',
                     'statistics_mc',
                     'T2KReWeight_flux'], COLORS['T2KReWeight_flux'])
FLUXSTAT['statistical (mc)'] = (['statistics_data',
                                 'statistics_mc'], COLORS['statistics_mc'])
FLUXSTAT['statistical (data)'] = (['statistics_data'], COLORS['statistics_data'])

T2KRW = ordereddict.OrderedDict()
T2KRW['flux'] = (['T2KReWeight_fsi',
                  'T2KReWeight_xs',
                  'T2KReWeight_flux'], COLORS['T2KReWeight_flux'])
T2KRW['cross-section'] = (['T2KReWeight_fsi',
                           'T2KReWeight_xs'], COLORS['T2KReWeight_xs'])
T2KRW['fsi'] = (['T2KReWeight_fsi'], COLORS['T2KReWeight_fsi'])

T2KFLUX = ordereddict.OrderedDict()
T2KFLUX['flux'] = (['T2KReWeight_flux'], COLORS['T2KReWeight_flux'])

T2KFSI = ordereddict.OrderedDict()
T2KFSI['fsi'] = (['T2KReWeight_fsi'], COLORS['T2KReWeight_fsi'])

T2KXS = ordereddict.OrderedDict()
T2KXS['cross-section'] = (['T2KReWeight_xs'], COLORS['T2KReWeight_xs'])

DETECTOR = ordereddict.OrderedDict()
DETECTOR['detector'] = (['highland_all_syst'], COLORS['detector'])

MASS = ordereddict.OrderedDict()
MASS['mass'] = (['Mass_all'], COLORS['Mass_all'])

OOFV = ordereddict.OrderedDict()
OOFV['oofv'] = (['highland_p0dccqeoofv_syst'], COLORS['detector'])

# for merged flux, xs, fsi throws
ALL = ordereddict.OrderedDict()
ALL['detector'] = (['statistics_data',
                    'statistics_mc',
                    'T2KReWeight_all',
                    'Mass_all',
                    'highland_all_syst'], COLORS['detector'])
ALL['mass'] = (['statistics_data',
                'statistics_mc',
                'T2KReWeight_all',
                'Mass_all'], COLORS['Mass_all'])
ALL['theory'] = (['statistics_data',
                  'statistics_mc',
                  'T2KReWeight_all'], COLORS['T2KReWeight_fsi'])
ALL['statistical (mc)'] = (['statistics_data',
                               'statistics_mc'], COLORS['statistics_mc'])
ALL['statistical (data)'] = (['statistics_data'], COLORS['statistics_data'])

# for ungrouped t2krw dials
UNGROUPEDXS = ordereddict.OrderedDict()
UNGROUPEDXS['MaQE'] = (['T2KReWeight_NIWG2014a_pF_O16',
                        'T2KReWeight_NIWG2014a_Eb_O16',
                        'T2KReWeight_NIWGMEC_Norm_O16',
                        'T2KReWeight_NXSec_CA5RES',
                        'T2KReWeight_NXSec_MaNFFRES',
                        'T2KReWeight_NXSec_MaCCQE'], COLORS['T2KReWeight_xs'])
UNGROUPEDXS['MaNFFRES'] = (['T2KReWeight_NIWG2014a_pF_O16',
                        'T2KReWeight_NIWG2014a_Eb_O16',
                        'T2KReWeight_NIWGMEC_Norm_O16',
                        'T2KReWeight_NXSec_CA5RES',
                        'T2KReWeight_NXSec_MaNFFRES'], COLORS['T2KReWeight_fsi'])
UNGROUPEDXS['CA5RES'] = (['T2KReWeight_NIWG2014a_pF_O16',
                        'T2KReWeight_NIWG2014a_Eb_O16',
                        'T2KReWeight_NIWGMEC_Norm_O16',
                        'T2KReWeight_NXSec_CA5RES'], COLORS['T2KReWeight_flux'])
UNGROUPEDXS['MEC_O16'] = (['T2KReWeight_NIWG2014a_pF_O16',
                        'T2KReWeight_NIWG2014a_Eb_O16',
                        'T2KReWeight_NIWGMEC_Norm_O16'], COLORS['detector'])
UNGROUPEDXS['Eb_O16'] = (['T2KReWeight_NIWG2014a_pF_O16',
                        'T2KReWeight_NIWG2014a_Eb_O16'], COLORS['statistics_mc'])
UNGROUPEDXS['pF_O16'] = (['T2KReWeight_NIWG2014a_pF_O16'], COLORS['statistics_data'])

MAQE = ordereddict.OrderedDict()
MAQE['MaQE'] = (['T2KReWeight_NXSec_MaCCQE'], COLORS['T2KReWeight_xs'])
EBO = ordereddict.OrderedDict()
EBO['Eb_O16'] = (['T2KReWeight_NIWG2014a_Eb_O16'], COLORS['statistics_mc'])
PFO = ordereddict.OrderedDict()
PFO['pF_O16'] = (['T2KReWeight_NIWG2014a_pF_O16'], COLORS['statistics_data'])
MECO = ordereddict.OrderedDict()
MECO['MEC_O16'] = (['T2KReWeight_NIWGMEC_Norm_O16'], COLORS['detector'])
##################### End #######################

#################### Volumes ####################
DEFAULT_VOLUME = 'Entire P0D-FV'
XLAYER_VOLUME = 'X-layer Only'
YLAYER_VOLUME = 'Y-layer Only'
BAGXLAYER_VOLUME = 'No Y-layers'
##################### End #######################

################### Reweight ####################
BANFF_FILE = os.path.join(directories.IRODS,
                          'asg/banff/postfit/postfit_banff_2015_data_20150417_allparams.root')
MASS_FILE = os.path.join(os.getenv('P0DNUMUCCANALYSISROOT'), 'data', 'P0DMassWeights.dat')
EMC_FILE = os.path.join(directories.WORK314, 'myAnalyses', 'P0DCC0Pi', 'Systematics', 'emc.root')
##################### End #######################

################### xsCache #####################
XSCACHE_FILE = os.path.join(directories.XSCACHE_STORE, 'xscache.root')
##################### End #######################

############# flux test files ###################
# these files were generated with 5000 throws of the t2kreweight
# parameters using a small highland2 microtree. this allows us to get
# a large number of throws to calculate the flux uncertainties.
JEREMY_TEST = '/usr/users/jplopez/data41/T2K/Work/P0DCC0Pi/FluxTest/'
MCP_TEST = os.path.join(JEREMY_TEST, 'microtree.root')
RDP_TEST = os.path.join(JEREMY_TEST, 'microtree.root')
WEIGHTS_TEST = os.path.join(JEREMY_TEST, 'out2', 'weights_unified.root')
RWFILE_DEFAULT_TEST = os.path.join(JEREMY_TEST, 'out2', 'rwfile_default.root')
RWFILE_TRUTH_TEST = os.path.join(JEREMY_TEST, 'out2', 'rwfile_truth.root')
##################### End #######################

### DEBUG ###
if DEBUG:
    # Some alternative selections
    SELECT_XLAYER = SELECTIONPAIR(SELECT.reco+'&&'+getLayerCutStr('selmu_pos[2]', 'x'),
                                  SELECT.true+'&&'+getLayerCutStr('truemu_pos[2]', 'x'))
    SELECT_YLAYER = SELECTIONPAIR(SELECT.reco+'&&'+getLayerCutStr('selmu_pos[2]', 'y'),
                                  SELECT.true+'&&'+getLayerCutStr('truemu_pos[2]', 'y'))
    SELECT_BAGX = SELECTIONPAIR(SELECT.reco+'&&'+getLayerCutStr('selmu_pos[2]', 'bagx'),
                                SELECT.true+'&&'+getLayerCutStr('truemu_pos[2]', 'bagx'))

    # testing 100% pure sideband
    PURE_SIDEBANDS = [SELECTIONPAIR('Sum$(accum_level[0][1]) > 3 && np0d == 2 && topology==1', 'topology == 1'),
                      SELECTIONPAIR('Sum$(accum_level[0][1]) > 3 && np0d > 2 && topology==2', 'topology == 2')]

    MICHEL1_SIDEBANDS = [SELECTIONPAIR('Sum$(accum_level[0][1]) > 3 && np0d == 2 && nmichel==1', 'topology == 1'),
                         SELECTIONPAIR('Sum$(accum_level[0][1]) > 3 && np0d > 2', 'topology == 2')]

    # testing sideband with additiona dphiT cut
    DPHIT_SIDEBANDS = [SELECTIONPAIR(
        'Sum$(accum_level[0][1]) > 3 && np0d == 2 && sel_dphi*180/TMath::Pi()>60', 'topology == 1')]

    # debugging source error dependece on the order of sources, fixed due to bincache
    # DETECTOR['detector'] = (['T2KReWeight_fsi', 'highland_all_syst'], COLORS['detector'])
