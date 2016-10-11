#!/usr/bin/env python

import sys
from tianlu.production.lib import MCConfigClass, MCJobClass
from tianlu.production.lib.utils import parseMissingFileList

def main(): 
    runInfo = MCJobClass.MCRunInfoClass()
    runCfgInfo = MCConfigClass.MCConfigClass(runInfo)
    print runInfo.nd280Version
    if len(sys.argv) == 1:
        runCfgInfo.Usage()
        sys.exit()
    if len(sys.argv) == 2:
        arg = sys.argv[1]
        if arg.title() == 'Help':
            Usage()
            sys.exit()
        if arg.title() == 'hardcoded':
            for run in range(0,8):
                for subrun in range(0, 100):
                    runCfgInfo.writeCfgFile(run, subrun)
        if '.pkl' in sys.argv[1]:
            missingDict = parseMissingFileList(sys.argv[1])
            for run in missingDict.keys():
                for subrun in missingDict[run]:
                    runCfgInfo.writeCfgFile(run%1000,subrun)
    if len(sys.argv) == 3:
        runs = runCfgInfo.parseRunSubrunArgs(sys.argv[1])
        subs = runCfgInfo.parseRunSubrunArgs(sys.argv[2])
        for runN in runs:
            for subN in subs:
                runCfgInfo.writeCfgFile(runN, subN)
    if len(runInfo.MissingConfigs) > 0:
        print 'Missing input file summary:'
        print runInfo.MissingConfigs

if __name__ == '__main__':
    main()
