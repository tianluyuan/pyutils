"""
Binning definitions for variables variables
"""
import array
import math
from tianlu.analysis.P0DCutUtils import getZBinLowX
from tianlu import rootutils as rtutils


#################### Binning ####################
MOM = (array.array("d", [0, 0.4, 0.5, 0.7, 0.9, 2.5, 5]),)
# Taken from Melody's thesis [units GeV]
MOM_MELODY = (array.array("d", [0, 0.4, 0.5, 0.7, 0.9, 2.5, 5]),)
# From Andy F. CC0pi TN215
MOM_ANDY = (array.array('d', [0, 0.2, 0.35, 0.5, 0.65, 0.8, 0.95, 1.1, 1.25,
                        1.5, 2, 3, 5]),)
MOM_FINE = (array.array('d', [i*0.1 for i in range(51)]),)
# Transverse momentum binning
MOMT = (array.array('d', [i*0.3 for i in range(8)]),)
# Q2
Q2 = (array.array('d', [i*0.09 for i in range(16)]),)
# Bjorken x
X = (array.array('d', [i*0.02 for i in range(51)]),)
# CosTheta
CTH = (array.array('d', [-1, 0, .6, .7, .8, .85, .9, .925, .975, 1]),)
CTH_PI2 = (array.array('d', [0, .6, .7, .8, .85, .9, .925, .975, 1]),)
CTH_PI2_SARA = (array.array('d', [0, .6, .7, .8, .85, .9, .94, .98, 1]),)
# Theta
TH = (array.array('d', [0, 5, 15, 30, 45, 90, 180]),)
TH_PI2 = (array.array('d', [0, 5, 15, 30, 45, 90]),)
TH_FINE = (array.array('d', [i for i in range(181)]),)
TH_INSET_FINE = (array.array('d', [i/60. for i in range(181)]),)
TH_LOGX = (array.array('d', [i*180/math.pi for i in rtutils.bin_log(10,
                                                                    0,
                                                                    math.pi,
                                                                    force_min=True)]),)
TH_LOGX_INSET = (TH_LOGX[0][0:5],)
TH_LOGX_NEG = (array.array('d', [i*-180/math.pi
                                for i in reversed(rtutils.bin_log(10, 0, math.pi))]),)
TH_LOGX_2PI = (TH_LOGX_NEG[0]+TH_LOGX[0],)
# categ binning
CATEG = (array.array('d', [i for i in range(11)]),)
# nparticles binning
NPARTICLES = (array.array('d', [i for i in range(11)]),)
NPARTICLES_EXTEND = (array.array('d', [i for i in range(51)]),)
# P0D z-binning
POSZ = (array.array("d", getZBinLowX()),)
POSZ_WTRBGS = (array.array("d", getZBinLowX(water_bags=True)),)
# email from Paul Rojas
POSZ_CSU = (array.array('d', [-3500, -3285, -3267, -3244, -3225, -3200,
                             -3183, -3156, -3139, -3113, -3097, -3070, -3053, -3027, -3005, -2970,
                             -2940, -2903, -2870, -2836, -2800, -2767, -2735, -2699, -2666, -2632,
                             -2600, -2564, -2530, -2497, -2463, -2429, -2395, -2360, -2329, -2292,
                             -2261, -2225, -2191, -2157, -2124, -2079, -2049, -2012, -1986, -1943,
                             -1912, -1875, -1849, -1808, -1779, -1740, -1710, -1673, -1640, -1605,
                             -1575, -1537, -1505, -1468, -1437, -1400, -1369, -1333, -1302, -1266,
                             -1240, -1207, -1190, -1165, -1146, -1122, -1105, -1078, -1059, -1034,
                             -1017, -987, -963, -947, -900]),)

# 2d
MOM_CTH_2D = (MOM[0], CTH[0])
MOM_CTH_PI2_2D = (MOM[0], CTH_PI2[0])
MOM_TH_PI2_2D = (MOM[0], TH_PI2[0])
MOM_TH_2D = (MOM[0], TH[0])
MOM_TH_LOGX_2D = (MOM[0], TH_LOGX[0])

MOM_MELODY_CTH_2D = (MOM_MELODY[0], CTH[0])
MOM_MELODY_CTH_PI2_2D = (MOM_MELODY[0], CTH_PI2[0])
MOM_MELODY_TH_2D = (MOM_MELODY[0], TH[0])
MOM_MELODY_TH_PI2_2D = (MOM_MELODY[0], TH_PI2[0])
MOM_MELODY_TH_LOGX_2D = (MOM_MELODY[0], TH_LOGX[0])

Q2_TH_LOGX_2D = (Q2[0], TH_LOGX[0])
MOMT_TH_LOGX_2D = (MOMT[0], TH_LOGX[0])
TOPOLOGY_REACTION_2D = (CATEG[0], CATEG[0])
##################### End #######################
