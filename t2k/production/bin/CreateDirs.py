#!/usr/bin/env python
#
#############################################################################
#
# Title: ND280 MCP directory structure initialization script
#
# Usage: This script should be called from the base production directory, i.e.
#                /nfs/wafflesciences/t2k/production005/
#
#                Ensure you have the desired generators, baselines, volumes and 
#                MC types set in the cardname. The cardname should have the following
#                format:
'''
production:         5
respin:                 F
nuMCType:               genie
verify:                 True
version:                v10r11p21
'''
#                the whole directory structure as defined here:
#                http://www.t2k.org/nd280/datacomp/howtoaccessdata/directorystructure 
#                and here:
#                http://www.t2k.org/nd280/datacomp/mcproductionruns/production004.
#                Note, the actual sample configurations are generated using 
#                conditionals within the loop below.
#
#                Finally, this script will link the configuration file creation scripts
#                and the job submission scripts into each sample directory, and you
#                should call them from there
#
#                This script was lovingly stolen and modified from Patrick DePerio's 
#                mk_dir.sh SciNet script
#############################################################################

import os
import sys

def fillFromCard(cardname):
    input = open(cardname, 'r')
    lines = []
    dirDict = {}
    for line in input:
        dirDict[line.split(':')[0]] = line.split(':')[1].strip()
    
    keyCheck = ['production', 'respin', 'nuMCType', 'verify', 'version']
    dirKeys = dirDict.keys()
    for key in keyCheck:
        if key not in dirKeys:
            print 'Check', cardname, 'for proper format, could not find \'%s\' key'%key
            sys.exit()
    return dirDict

def checkCombination(mcType, baseline, vol, beamType):
    if mcType == 'anti-genie':
        return (baseline == '2010-11-water' and
                ((vol == 'magnet' and
                 beamType in ['run5']) or
                (vol == 'basket' and
                 beamType in ['tpcgas'])))
    elif mcType == 'genie':
        if vol == 'magnet':
            if baseline == '2010-02-water':
                if beamType != 'run1':
                    return False # only run1 for 2010-02-water
                return True # 2010-02-water run1 magnet
            if baseline == '2010-11-air':
                if (beamType not in ['run2', 'run3', 'run4']):
                    return False # beamb for Run2, beamc for Run3/4
                return True # 2010-11-air has run2,3,4 magnet
            if baseline == '2010-11-water':
                if (beamType not in ['run2', 'run4', 'run5']):
                    return False # beamb for Run2, beamc for Run3/4, beamd for Run5
                return True # 2010-11-water run2,4,5 magnet

        if vol == 'basket':
            if '2010-11' not in baseline:
                return False # only 2010-11 for basket
            if 'run' in beamType:
                return False # no full spill 'runN'
            return True # 2010-11 water/air basket {beam, nue, pistuff)
    else:
        print 'Undefined nuMCType', mcType
        print 'Will not make dirs!'
        return False

def makePath(aList):
    path = ''
    for item in aList:
        path += str(item)+'/'
    path = path[0:-1] # strip out the last /
    return path

def main():
    createDirs = False
    if len(sys.argv) == 2:
        cardName = sys.argv[1]
    else:
        print 'Need to provide the cardname'
        sys.exit()
    inp = raw_input('Create Directories for real?:')
    if inp.title() == 'Y':
        createDirs = True
    else:
        createDirs = False

    dirDict = fillFromCard(cardName)
    #god help you if you don't run this not from the    production directory
    prodDir = os.environ['PWD'] 
    print 'Production Directory is', prodDir

    copiedFiles = ['runInfo.card']
    respin = dirDict['respin']
    simName = dirDict['nuMCType']
    version = dirDict['version']

    baselines = ['2010-02-water', '2010-11-water', '2010-11-air']
    volumes = ['magnet', 'basket']
    beamTypes = ['run1',
                 'run2',
                 'run3',
                 'run4',
                 'run5',
                 'beam',
                 'ccpiplus',
                 'ccpizero',
                 'ncpiplus',
                 'ncpizero',
                 'nue',
                 'tpcgas']

    if dirDict['verify'] == 'True':
        simName = 'verify/'+version+'/'+simName
        volumes = ['magnet', 'basket']
        # Don't need all beamTypes for verification
        beamTypes = ['beam', 'run2', 'run3', 'run4']

    for baseline in baselines:
        for volume in volumes:
            for beamType in beamTypes:
                if checkCombination(simName, baseline, volume, beamType):
                    print '\nCreating file for:\n', respin, simName, baseline, volume,
                    print beamType
                    # a good combination, so let's make a directory!
                    subdirs=['anal','cali','dats','elmc','g4mc','gnmc','numc','nucp','reco','staging','logf']
                    dirs = [makePath([prodDir,respin])]
                    for subdir in subdirs:
                        dirs.append(makePath([prodDir, respin, 'mcp', simName, baseline,
                                                volume, beamType, subdir]))
                    for step in ['nd280', 'nucp', 'numc']:
                        for substep in ['logs', 'cfg', 'scripts']:
                            subdir = step+'/'+substep
                            path = makePath([prodDir,respin,'mcp',simName,baseline,volume, 
                                                            beamType, 'staging', subdir])
                            dirs.append(path)
                    for dir in dirs:
                        cmd = 'mkdir -m 775 -p '+dir
                        print cmd
                        if createDirs:
                            os.system(cmd)
                        if dir[-7:] == 'staging':
                            cmds = []
                            for file in copiedFiles:
                                cmds.append('cp '+prodDir+'/'+file+' '+dir)
                            for cmd in cmds:
                                print cmd
                                if createDirs:  
                                    os.system(cmd)

if __name__ == '__main__':
    main()
