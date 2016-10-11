#!/usr/bin/env python

import os
import sys
import tianlu.production.lib.utils as produtils
from optparse import OptionParser

# add a cmdline parser to take optional directory specifications
usage = 'usage: %prog [options] path_to_rundir run subruns[or missing.pkl]'
parser = OptionParser(usage)
parser.add_option('-a', '--all',
                  dest='purge_all',
                  action='store_true',
                  help='Look for files_to_purge in all MC directories')
parser.add_option('-g', '--genie',
                  dest='purge_genie',
                  action='store_true',
                  help='Look for files_to_purge in genie MC directories')
parser.add_option('-c', '--nucp',
                  dest='purge_nucp',
                  action='store_true',
                  help='Look for files_to_purge in nucp directory')
parser.add_option('-k', '--key',
                  dest='specifier',
                  action='store',
                  default='missing',
                  type='string',
                  help='Specify key to purge run-subrun from pkl_dict[key].')

(opt, args) = parser.parse_args()
base_run_dir = args[0]

if args[1].endswith('.pkl'):
    try:
        missingDict = produtils.parseMissingFileList(args[1], opt.specifier)
    except KeyError:
        print 'no such key,', opt.specifier, 'in missing.pkl!'
        sys.exit(1)
else:
    run = int(args[1])
    subruns = []
    for arg in args[2:]:
        subruns.append(int(arg))
    missingDict = {run:subruns}

# set which directories to purge
dirs_genie = ['gnmc', 'numc']
dirs_nucp = ['nucp']
dirs_nd280 = ['anal', 'cali', 'elmc', 'g4mc', 'reco']
if opt.purge_genie:
    dirs = dirs_genie
elif opt.purge_nucp:
    dirs = dirs_nucp
elif opt.purge_all:
    dirs = dirs_genie + dirs_nd280 + dirs_nucp
else:
    dirs = dirs_nd280

toPurge = []
for a_dir in dirs:
    full_dir = os.path.join(base_run_dir, a_dir)
    for a_file in os.walk(full_dir).next()[-1]:
        run, subrun = produtils.get_run_and_subrun(a_file)
        if missingDict.has_key(run):
            if subrun in missingDict[run]:
                toPurge.append(os.path.join(full_dir, a_file))
        elif missingDict.has_key(run%1000):
            if subrun in missingDict[run%1000]:
                toPurge.append(os.path.join(full_dir, a_file))

print 'These are the files that will be deleted!'
for file in toPurge:
    print file

doIt = raw_input('Do you want to purge the above files??? [y/n]')
if doIt == 'y' or doIt == 'Y':
    for file in toPurge:
        cmd = 'rm '+file
        print cmd
        os.system(cmd)


