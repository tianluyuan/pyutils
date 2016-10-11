""" Module containing settings for data directories
"""
# work directories
WORK304 = '/axp/lnxdata304/wrk2/ltm/tianlu/'
WORK314 = '/nfs/data44/t2k/tianlu/'
SOFTWARE = WORK314+'software/'
HADOOP_T2K = '/mnt/hadoop/store/user/tianlu/t2k/'

# data store directories for downloaded data
HADOOP_MARINO_GRID = '/mnt/hadoop/store/user/marino-grid/'
HADOOP_JPLOPEZ_GRID = '/mnt/hadoop/store/user/jplopez-grid/'
DATA41_T2K = '/nfs/data41/t2k/'
DATA41_JPLOPEZ = '/nfs/data41/t2k/jplopez'
IRODS = DATA41_T2K+'irods/'

# post-runsort data directories
NFS_BASE = WORK314+'database/'
HADOOP_BASE = HADOOP_T2K+'data/'

# paths for working with grid
LFN_BASE = '/grid/t2k.org/'
LFN_ND280 = LFN_BASE+'nd280/'
CU_SRM01_SERVER = 'srm://hepse01.colorado.edu:8443/srm/v2/server'
CU_SRM02_SERVER = 'srm://hepse02.colorado.edu:8443/srm/v2/server'
TRIUMF_SRM = 'srm://t2ksrm.nd280.org'

# production directories
HADOOP_PROD = HADOOP_T2K+'GNProd/'
NFS_DATA44_PROD = '/nfs/data44/t2k/GNProd/'
NFS_DATA42_PROD = '/nfs/data42/t2k/GNProd/'
ALL_PRODUCTION_DIRS = [HADOOP_PROD,
                       NFS_DATA44_PROD,
                       NFS_DATA42_PROD]

# analysis output directories
STORE = WORK314+'store/'
HIGHLAND_STORE = STORE+'highland/analyses/'
HIGHLAND2_STORE = STORE+'highland2/analyses/'
RUNANALYSES_STORE = STORE+'other/runanalyses/'
TIANND280_STORE = STORE+'other/tiannd280/'
WEIGHTS_HIGHLAND_STORE = STORE+'xstool/weights_highland/'
WEIGHTS_HIGHLAND2_STORE = STORE+'xstool/weights_highland2/'
WEIGHTS_POSTFIT_STORE = STORE+'xstool/weights_postfit/'
WEIGHTS_UNGROUPED_STORE = STORE+'xstool/weights_ungrouped/'
XSCACHE_STORE = STORE+'xstool/cache/'
GENERATORS_STORE = STORE+'other/generators/'
TRASH_STORE = WORK314+'.trash/'

#### Data locations for downloading ####
# this may change if we want to download files to different disk
# Primary locations are HADOOP_MARINO_GRID, HADOOP_JPLOPEZ_GRID, or DATA41_T2K
DATASTORE = HADOOP_MARINO_GRID
