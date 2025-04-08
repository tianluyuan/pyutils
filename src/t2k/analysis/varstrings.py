"""
List of variables that we may want to plot for analyses.
Variables are formed using TTree branch names.

Sum$ is to force ROOT to evaluate the formula. See xsTutorial.dox
"""
from tianlu.analysis.constants import *


################### highland ######################
TRUE_NU_E = 'nu_trueE/1e3'  # GeV
TRUE_MU_MOM = 'truemu_mom/1e3'
TRUE_MU_E = 'sqrt(({mu_p})**2+({mu_m})**2)'.format(mu_p=TRUE_MU_MOM, mu_m=MU_MASS)
TRUE_MU_POSZ = 'truemu_pos[2]'
TRUE_MU_CTH = 'truemu_costheta'
TRUE_MU_TH = 'TMath::ACos(truemu_costheta)*180/TMath::Pi()'
TRUE_Q2 = '2*{nu_e}*({mu_e}-{mu_p}*{mu_cth})-{mu_m}**2'.format(mu_m=MU_MASS,
                                                               nu_e=TRUE_NU_E,
                                                               mu_e=TRUE_MU_E,
                                                               mu_p=TRUE_MU_MOM,
                                                               mu_cth=TRUE_MU_CTH) # approximation
TRUE_X = '{Q2}/(2*{n_m}*({nu_e}-{mu_e}))'.format(Q2=TRUE_Q2,
                                                 n_m=(P_MASS+N_MASS)/2,
                                                 nu_e=TRUE_NU_E,
                                                 mu_e=TRUE_MU_E)

TRUE_DPHI = 'true_dphi*180/TMath::Pi()'
TRUE_MU_MOMT = 'truemu_momt/1e3'
TRUE_P_MOMT = 'truep_momt/1e3'
SEL_DPHI = 'sel_dphi*180/TMath::Pi()'
SEL_MU_MOM = 'Sum$(selmu_mom[0])/1e3'
SEL_MU_MOMT = 'selmu_momt/1e3'
SEL_MU_CTH = 'TMath::Cos(Sum$(selmu_theta[0]))'
SEL_MU_TH = 'Sum$(selmu_theta[0])*180/TMath::Pi()'
SEL_MU_POSZ = 'selmu_pos[2]'
DELTA_DPHI = '{0}-{1}'.format(SEL_DPHI, TRUE_DPHI)

# particle counters
NPROTON = 'nproton'
NNEUTRON = 'nneutron'
NMESON = 'nmeson'
NPHOTON = 'nphoton'
NMICHEL = 'nmichel'

# 2d
TRUE_MU_MOM_CTH_2D = '{1}:{0}'.format(TRUE_MU_MOM, TRUE_MU_CTH)
TRUE_MU_MOM_TH_2D = '{1}:{0}'.format(TRUE_MU_MOM, TRUE_MU_TH)
TRUE_TOPOLOGY_REACTION_2D = '{1}:{0}'.format('topology', 'reaction')
TRUE_NU_E_DPHI_2D = '{1}:{0}'.format(TRUE_NU_E, TRUE_DPHI)
TRUE_Q2_DPHI_2D = '{1}:{0}'.format(TRUE_Q2, TRUE_DPHI)
TRUE_MU_MOMT_DPHI_2D = '{1}:{0}'.format(TRUE_MU_MOMT, TRUE_DPHI)
SEL_MU_MOM_CTH_2D = '{1}:{0}'.format(SEL_MU_MOM, SEL_MU_CTH)
SEL_MU_MOM_TH_2D = '{1}:{0}'.format(SEL_MU_MOM, SEL_MU_TH)
SEL_MU_MOM_DPHI_2D = '{1}:{0}'.format(SEL_MU_MOM, SEL_DPHI)
##################### End #######################
