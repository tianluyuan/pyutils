#!/usr/bin/env python

import sys
import subprocess


def get_volume(file, baseFile):
    # naming convention different for 13a
    if '13a_nom_250ka' in file:
        volume = 'nd5'
    elif '13a_nom_m250ka' in file:
        volume = 'nd5m'
    elif '13a_nom_ND6_250ka' in file:
        volume = 'nd6'
    elif '13a_nom_ND6_m250ka' in file:
        volume = 'nd6m'
    else:
        # otherwise we're not using 13a flux
        volume = baseFile[baseFile.find('.') + 1:baseFile.find('_')]

    return volume


def checkNameAndVolumeAndSort(fileList):
    first = True
    for file in fileList:
        if file.find('/') != -1:
            # just grab the bare file name (sans path)
            baseFile = file[file.rfind('/') + 1:]
        else:
            baseFile = file
        if 'fluka_13a' not in baseFile and '_flukain' not in baseFile:
            print basefile, 'doesn\'t have \'fluka\''
            return -1
        if baseFile[-5:] != '.root':
            print baseFile, 'is not a ROOT file'
            return -1
        # grab 'nd9' from 'nu.nd9_flukain.xxx.root
        if first:
            volume = get_volume(file, baseFile)
            first = False

        currentVolume = get_volume(file, baseFile)
        if volume != currentVolume:
            # then the user has given a list with a mix of geometries
            print 'Found baseFiles with different volumes:', volume, currentVolume
            print 'Please provide a list with a single volume'
            return -1
    # sort the file list in order of file number
    # fileList.sort(key=lambda n: int(n[-10:-5][n[-10:-5].rfind('.')+1:]))
    fileList.sort()

    return volume


def main():
    usage = '\n USAGE: ' + sys.argv[0] + ' <list of jNuBeam flux files>\n'
    usage += 'i.e. ' + sys.argv[0] + ' path/to/files/nu.nd5_flukain.*.root\n'
    if len(sys.argv) == 1:
        print usage
        return 0

    filesToHadd = sys.argv[1:]
    print 'files to hadd:', filesToHadd
    volume = checkNameAndVolumeAndSort(filesToHadd)
    cmd = 'hadd -f nu.' + volume + '_flukain.0-' + \
        str(len(filesToHadd) - 1) + '.root ' + ' '.join(filesToHadd)
    # for file in filesToHadd:
    #   cmd += ' ' + file
    # print cmd
    history = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE).communicate()
    if history[1] == '':
        print history[0]
    else:
        print 'Log:'
        print history[0]
        print 'Errors:'
        print history[1]
if __name__ == '__main__':
    main()
