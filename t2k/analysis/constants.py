"""
Provides constants to be used in analysis.  Includes both physics
and detector constants.
"""


__author__ = 'Tianlu Yuan'
__email__ = 'tianlu.yuan [at] colorado.edu'


### Physics###
MU_MASS = 0.1057
P_MASS = 0.9383
N_MASS = 0.9396

### Detector ###
# new P0D Fiducial from Karin P0DGeom (Erez uses)
XCUT_FIDUCIAL_MIN = -836
XCUT_FIDUCIAL_MAX = 764
YCUT_FIDUCIAL_MIN = -871
YCUT_FIDUCIAL_MAX = 869
ZCUT_FIDUCIAL_MIN = -2969
ZCUT_FIDUCIAL_MAX = -1264
# ecal p0dule z
Z_EC_WIDTH = 44
# water p0dule z
Z_WT_WIDTH = 68
# water bag z
Z_BG_WIDTH = 29
# p0dule width
Z_PD_WIDTH = Z_WT_WIDTH-Z_BG_WIDTH
# USWT has a water cover of 6.35 mm thick to support last bag
Z_WC_WIDTH = 7
# defining z positions
# TN73 Karin's FV definition goes to half-way in the first upstream
# water p0dule.  We define wt_upstream to be midway between USECAL and
# USWT.  The entire z-binning is based off this coordinate
WT_UPSTREAM = -2969-Z_PD_WIDTH/2
# number of p0dules upstream/downstream of fv_upstream to bin
N_UPSTREAM_P0DULES = 8
N_UPSTREAM_P0DULES_sand = 11
N_DOWNSTREAM_P0DULES = 5
# number of water bags in each USWT and CWT.  used when want to bin water
N_BGS_USWT = 13
N_BGS_CWT = 12
# Cutoff for bunch-0 or before-beam-arrives bunch [ns]
T_BUNCH0_CUTOFF = 2500
# threshold for matching p0d-tpc tracks
TPC1FMOM_MIN = 0
