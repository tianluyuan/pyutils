#!/usr/bin/env python

import sys
import pickle
import re
import math
import itertools
import glob
import os
from collections import defaultdict
from optparse import OptionParser
from tianlu.production.lib.utils import get_run_and_subrun


# add a cmdline parser to take optional n-subrun values
usage = 'usage: %prog [options] path_to_dir'
parser = OptionParser(usage)
parser.add_option('-s',
                  type='int',
                  dest='n_subrun',
                  help='Set the number of subruns that should exist.')
parser.add_option('-o',
                  default='Missing.pkl',
                  type='string',
                  dest='output',
                  help='Set the output pkl file containing missing files.')

(opt, args) = parser.parse_args()


def check_handle(handle):
    """ regular expressions check to make sure filename qualifies for parsing
    """
    rex = r'oa_.+_[0-9]{8}-[0-9]{4}_.+\.root'
    return re.match(rex, handle)


filedir = args[0]
outfilename = opt.output

runs = defaultdict(list)
missing = {}
extra = {}
anoms = defaultdict(list)
anoms_list = []
duplicates = []
fileSizes  = []
checkSizes = True
print filedir
lsFileList = glob.glob(os.path.join(filedir, '*.root'))

for a_file in lsFileList:
    filename = os.path.basename(a_file)
    if not check_handle(filename):
        print filename, 'of incorrect run-subrun format.  Skipped!'
        continue
    size = os.path.getsize(a_file)
    fileSizes.append((size, filename))
    if 'root' not in filename:
        print 'wtf?', filename
        sys.exit()
    run, subrun = get_run_and_subrun(filename)
    if subrun in runs[run]:
        duplicates.append((run,subrun))
    runs[run].append(subrun)

if len(runs) == 0:
    print 'No run-subrun files in dir.'
    sys.exit(1)

nExpected = 0
nMissing = 0
nExtra = 0
# Set nSubRuns to default to the max_subrun+1 in the first five runs
# this should be enough to guarantee we have the correct nSubruns
if not opt.n_subrun:
    nSubRuns = max(itertools.chain(*runs.values()[0:4]))+1
else:
    nSubRuns = opt.n_subrun

# use max_run for looping in case an entire run is missing
min_run = min(runs)
max_run = max(runs)
print 'Checking run ranges', min_run, 'to', max_run
for run in range(min_run, max_run+1):
    thisRun = runs[run]
    # using set difference, which will lose duplicates when diffing but
    # this is what we want because we don't want to purge away dupes
    amiss = set(range(nSubRuns))-set(thisRun)
    aextra = set(thisRun)-set(range(nSubRuns))

    if len(amiss) != 0:
        print '--------------------------------------------------------'
        print 'something is amiss with this run', run
        print 'Number of elements', len(thisRun)
        print thisRun
        missing[run] = amiss

    if len(aextra) != 0:
        print '--------------------------------------------------------'
        print 'something is extra with this run', run
        print 'Number of elements', len(thisRun)
        print thisRun
        extra[run] = aextra

    nMissing += len(amiss)
    nExtra += len(aextra)
    nExpected += nSubRuns

print '%(miss)i missing out of %(exp)i expected, or %(div).2f%% completion' %\
                {"miss":nMissing, "exp":nExpected, "div":100-100*float(nMissing)/nExpected}
for run in missing.keys():
    print run, missing[run]

print '%(extra)i extra out of %(exp)i expected, or %(div).2f%% over' %\
                {"extra":nExtra, "exp":nExpected, "div":100*float(nExtra)/nExpected}
for run in extra.keys():
    print run, extra[run]

if len(duplicates) > 0:
    print 'Found Duplicates:'
    for dupe in duplicates:
        print dupe
    print 'Please please please take care of duplicates manually before purging'

sizes, files = zip(*fileSizes)
avgSize = sum(sizes)/len(sizes)
cutoff = 0.15
for size, file in fileSizes:
    if math.fabs(size-avgSize)/float(avgSize) > cutoff: 
        anoms_list.append((file, size))
        print '  Size of file', file, 'is >',float(100*cutoff),'% different from avg file size!'
        print '      Avg:', avgSize, 'this file size:', size

        run, subrun = get_run_and_subrun(file)
        anoms[run].append(subrun)

# pickle missing, extra, and anomalous
outFile = open(outfilename, 'w')
pickle.dump({'missing':missing,
             'extra':extra,
             'anoms':anoms,
             'anoms_list':anoms_list}, outFile)
outFile.close()

