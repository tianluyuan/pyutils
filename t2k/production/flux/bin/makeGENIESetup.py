#! /usr/bin/env python
"""
create cfg files from template.cfg for rdp MIDAS processing
T.Y. 8/7/13
Modified for genie_flux_probs 11/8/13
Builds the cfg for runND280 and the script to submit to cluster
Usage: makeGENIESetup.py path/to/setup.py
where setup.py conatins the geometeries, flux_regions, and
master_volumes in lists.
"""
import os
import subprocess
from tianlu import utils as utl
import tianlu.production.flux as flux
from tianlu.production.lib import utils
import sys


__author__ = 'Tianlu Yuan'
__email__ = 'tianlu.yuan [at] colorado.edu'


cwd = utl.get_cwd()
template_dir = os.path.dirname(flux.__file__)
cfg_dir = os.path.join(cwd, 'cfg/')
sh_dir = os.path.join(cwd, 'scripts/')
evtrate_dir = os.path.join(cwd, 'genie_evtrate/')
log_dir = os.path.join(cwd, 'log/')
utl.make_dirs_if_needed(cfg_dir)
utl.make_dirs_if_needed(sh_dir)
utl.make_dirs_if_needed(evtrate_dir)
utl.make_dirs_if_needed(log_dir)
# This will be the output .cfg file
cfg_file_base = os.path.join(cfg_dir, 'genie_setup_flux_probs_{0}_{1}_{2}.cfg')
cfg_template = os.path.join(template_dir, 'templates', 'genie_setup_flux_probs.cfg')
sh_file_base = os.path.join(sh_dir, 'genie_setup_flux_probs_{0}_{1}_{2}.sh')
sh_template = os.path.join(template_dir, 'templates', 'template_script.sh')

# dynamic import of command line arg
sys.path.append(os.path.realpath(os.path.dirname(sys.argv[1])))
setup = __import__(os.path.splitext(sys.argv[1])[0], globals(), locals())

def create_dict(geometry, flux_region, master_volume):
    cfg_dict = {}

    cfg_dict['%MASTER_VOLUME%'] = master_volume
    if 'nd6' in flux_region:
        cfg_dict['%FLUX_REGION%'] = 'Magnet'
    elif 'nd5' in flux_region:
        cfg_dict['%FLUX_REGION%'] = 'Basket'
    elif 'half' in flux_region:
        cfg_dict['%FLUX_REGION%'] = 'Basket'

    cfg_dict['%FLUX_FILE%'] = utl.glob_newest(cwd+'/nu.'+flux_region+'_flukain*')

    if 'water' in geometry:
        cfg_dict['%P0D_WATER%']='1'
    else:
        cfg_dict['%P0D_WATER%']='0'

    if geometry in utils.STANDALONE_GEOMETRIES:
        cfg_dict['%GEOMETRY%'] = 'standalone = '+geometry
    else:
        cfg_dict['%GEOMETRY%'] = 'baseline = '+geometry[0:7]

    cfg_dict['%N_FLUX_FILES%'] = str(int(cfg_dict['%FLUX_FILE%'].split('-')[-1].split('.root')[0])+1)
    cfg_dict['%XS_TABLE%'] = utl.glob_newest(cwd+'/gxspl-t2k-v*.xml')

    flux_probs_name = os.path.splitext(
        os.path.basename(cfg_dict['%FLUX_FILE%']))[0]
    flux_probs_name += '_genie_evtrate_'+geometry.lower()+'_'+master_volume.lower()+'.root'
    cfg_dict['%FLUX_PROBS_NAME%'] = flux_probs_name

    return cfg_dict

for geometry in setup.geometries:
    for flux_region in setup.flux_regions:
        for master_volume in setup.master_volumes:
            # make cfg files
            cfg_file_name = cfg_file_base.format(flux_region,
                                                 geometry,
                                                 master_volume.lower())
            cfg_dict = create_dict(geometry, flux_region, master_volume)
            utl.find_and_replace(cfg_template, cfg_file_name, **cfg_dict)

            # make sh scripts
            sh_file_name = sh_file_base.format(flux_region,
                                               geometry,
                                               master_volume.lower())
            sh_dict = {'%CFG_PATH%':cfg_file_name,
                       '%LOG%':log_dir,
                       '%EVTRATE%':evtrate_dir}
            utl.find_and_replace(sh_template, sh_file_name, **sh_dict)
            subprocess.call(['chmod', '755', sh_file_name])
