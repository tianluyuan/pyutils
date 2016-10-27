#!/usr/bin/env python
"""
copied from /nfs/data43/t2k/GNprod5/PurgeFromGrid.py by tianlu 6/17/14
"""
import subprocess
import os
import sys

pwd = os.environ['PWD']


if pwd[-7:] != 'staging':
  print 'You really should be running this from the \'staging\' directory'
  sys.exit()

run = int(sys.argv[1])
subrun = int(sys.argv[2])

dirs = ['anal', 'cali', 'elmc', 'g4mc', 'reco']
basePath = os.environ['PWD']
basePath = basePath[basePath.find('production'):]
basePath = basePath[0:basePath.find('staging')]
basePath = '/grid/t2k.org/nd280/'+basePath

toPurge = []

runHandle = '%(run)04i-%(subrun)04i'%{"run":run, "subrun":subrun}
for dir in dirs:
  path = basePath+dir
  lsCmd = subprocess.Popen('lfc-ls '+path, shell=True,
                              stdin =subprocess.PIPE, 
                              stdout = subprocess.PIPE).communicate()
  if lsCmd[1] != None:
    print 'uh oh', lsCmd[1]
    sys.exit()
  files = lsCmd[0].split('\n')
  for file in files:
    if runHandle in file:
      toPurge.append(path+'/'+file)

cmds = []
print 'Delete these files?'
for file in toPurge:
  cmds.append('lcg-del -a lfn:'+file)
  print cmds[-1]

toDelete = raw_input('[y/n] ' )
if toDelete == 'y' or toDelete == 'Y':
  for cmd in cmds:
    os.system(cmd)


