import ROOT
import time
import os
import sys
import numpy as np
from array import array
from thirdparty import autovivification
from tianlu import rootutils
from tianlu.analysis.constants import *


__author__ = 'Tianlu Yuan'
__email__ = 'tianlu.yuan [at] colorado.edu'


###############################################################
# Cut selection strings used in plotting hists and such
tpcCut = 'TPC1Used && FGDCosmic'
tpc1fidCut = ('TPC1FPos.X()>' + str(XCUT_FIDUCIAL_MIN) +
              ' && TPC1FPos.X()<' + str(XCUT_FIDUCIAL_MAX) +
              ' && TPC1FPos.Y()>' + str(YCUT_FIDUCIAL_MIN) +
              ' && TPC1FPos.Y()<' + str(YCUT_FIDUCIAL_MAX) +
              ' && TPC1FPos.Z()<-755')
p0dgoingCut = tpcCut + '&& TPC1FMom>250 &&' + tpc1fidCut
p0dgoingNoMomCut = tpcCut + '&&' + tpc1fidCut
p0dfoundCut = 'nP0Ds==1'
p0dtrackCut = p0dfoundCut + '&&P0DCone.Mag()==0'
p0dshowerCut = p0dfoundCut + '&&P0DCone.Mag()>0'
bkwrdscosCut = 'TPC1FDir.Y()>0&&TPC1FDir.CosTheta()>0.95'
# P0DStopping checks that the most upstream reconstructed P0DPosition is
# within a fiducial x-y plane and has its z within the P0D
p0dstoppingCut = ('((P0DBPos.Z()<P0DFPos.Z()'\
                  '&&P0DBPos.Z()>' + str(ZCUT_FIDUCIAL_MIN) +
                  '&&P0DBPos.X()>' + str(XCUT_FIDUCIAL_MIN) +
                  '&&P0DBPos.X()<' + str(XCUT_FIDUCIAL_MAX) +
                  '&&P0DBPos.Y()>' + str(YCUT_FIDUCIAL_MIN) +
                  '&&P0DBPos.Y()<' + str(YCUT_FIDUCIAL_MAX) +
                  ')||(P0DBPos.Z()>P0DFPos.Z()'\
                  '&&P0DFPos.Z()>' + str(ZCUT_FIDUCIAL_MIN) +
                  '&&P0DFPos.X()>' + str(XCUT_FIDUCIAL_MIN) +
                  '&&P0DFPos.X()<' + str(XCUT_FIDUCIAL_MAX) +
                  '&&P0DFPos.Y()>' + str(YCUT_FIDUCIAL_MIN) +
                  '&&P0DFPos.Y()<' + str(YCUT_FIDUCIAL_MAX) + '))')

# ##############
# Target PDG selection strings
atomicZ = rootutils.atomZ_sel_str('TargetPDG')
carbon = '({0} == 6)'.format(atomicZ)
copper = '({0} == 29)'.format(atomicZ)
zinc = '({0} == 30)'.format(atomicZ)
oxygen = '({0} == 8)'.format(atomicZ)
###############################################################

# ##############
# Color index scheme; for when ROOT colors suck
unityColorIdx = 38
goodColorIndices = [1, 2, 4, 797, 904, 8, 9, 882, 630, 598, 792]
#badColorIndices = [3, 5, 6, 7, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, unityColorIdx]

############
# Selections
############

def inP0DFiducial(position):
    """ Values taken from:
    ~/software/ND280/p0dAnalysis/macros/TianluP0DAnalysisMacros/P0DNuMuAnalysis/cutUtils.H
    """
    if (XCUT_FIDUCIAL_MIN < position.X() < XCUT_FIDUCIAL_MAX and
            YCUT_FIDUCIAL_MIN < position.Y() < YCUT_FIDUCIAL_MAX and
            ZCUT_FIDUCIAL_MIN < position.Z() < ZCUT_FIDUCIAL_MAX):
        return True
    else:
        return False


def inSandRegion(pid):
    # Restrict to ~first layer of P0D USECAL
    zCut_sand_max = -3260

    # The GlobalRecon Front Position is compromised by P0DTrack/P0DShower issue
    # It should not be used in our cuts when we're trying to study the Track/Shower
    # efficiencies.  Instead, we loop through the hits in the PID and check to see
    # if the first hit has frontZ position upstream of some zMax.  We sort these hits
    # by their time and see if the first registered hit is upstream of
    # zCut_sand_max.
    nHitsSaved = pid.NHitsSaved
    if nHitsSaved == 0:
        return False

    # hitTuple = (zpos, time)
    hitTuples = []
    for hitIdx in range(nHitsSaved):
        aHit = pid.HitsSaved.At(hitIdx)

        posT = aHit.Time
        posX = aHit.Position.X()
        posY = aHit.Position.Y()
        posZ = aHit.Position.Z()

        hitTuples.append((posZ, posX, posY, posT))

    # sort by first element in tuple
    hitTuples.sort()
    # print 'New Hits....'
    # print hitTuples

    # Get FIRST hit in PID and check if it's zPos is upstream of zCut_sand_max
    if (XCUT_FIDUCIAL_MIN < hitTuples[0][1] < XCUT_FIDUCIAL_MAX and
        YCUT_FIDUCIAL_MIN < hitTuples[0][2] < YCUT_FIDUCIAL_MAX and
        hitTuples[0][0] < zCut_sand_max):
        # print hitTuples
        return True
    else:
        return False


def reactionStrToInt(reactionStr):
    """ Groups GENIE reaction string into following categories

    Reactions can be divided into two basic categories:
    1) CC, 2) NC
    These can then be further divided into subcategories of
    interaction modes
    1) QES, 2) RES, 3)DIS, 4) COH
    And the target nucleon can be either:
    1) n, 2) p, 3) other (e.g. C12 in case of COH)
    We use a three digit integer to represent all the different
    possibilities, with the order being:
    first digit = target_nucleon
    second digit = CC_NC
    third digit = interaction_mode
    In the case of a string that doesn't fit any of the
    possibilities, a "0" is returned
    """
    reactionInt = 0
    if 'QES' in reactionStr:
        reactionInt += 1
    elif 'RES' in reactionStr:
        reactionInt += 2
    elif 'DIS' in reactionStr:
        reactionInt += 3
    elif 'COH' in reactionStr:
        reactionInt += 4

    if 'CC' in reactionStr:
        reactionInt += 10
    elif 'NC' in reactionStr:
        reactionInt += 20

    if '2112' in reactionStr:
        reactionInt += 100
    elif '2212' in reactionStr:
        reactionInt += 200
    else:
        reactionInt += 300

    return reactionInt


def reactionCodeToInt(reactionCode):
    """ Groups events by NEUT reaction code into the basic
    categories below

    Reactions can be divided into two basic categories:
    1) CC, 2) NC
    These can then be further divided into subcategories of
    interaction modes
    1) QES, 2) RES, 3)DIS, 4)Multi-Pi
    And the target nucleon can be either:
    1) n, 2) p
    We use a three digit integer to represent all the different
    possibilities, with the order being:
    first digit = target_nucleon
    second digit = CC_NC
    third digit = interaction_mode
    In the case of a string that doesn't fit any of the
    possibilities, a "0" is returned
    """
    reactionInt = 0
    absCode = np.abs(reactionCode)
    if absCode == 1 or absCode > 50:
        # QES
        reactionInt += 1
    elif absCode == 26 or absCode == 46:
        # DIS
        reactionInt += 3
    elif absCode == 21 or absCode == 41:
        # Multi-pi
        reactionInt += 4
    else:
        # RES
        reactionInt += 2

    if absCode < 30:
        # CC
        reactionInt += 10
    else:
        # NC
        reactionInt += 20

    # For our plotting purposes we don't really care whether target nucleon was p or n
    # NEUT treatment of target nucleon is kind of a mess.  I'm being hacky here by
    # just treating all nucleon as a n.
    reactionInt += 100

    return reactionInt


def getRunNum(IsMC, RunID):
    if not IsMC:
        if RunID < 6462:
            return 1
        elif RunID < 8115:
            return 2
        elif RunID < 8983:
            return 3
        else:
            return 4
    else:
        return int(str(RunID)[2])


def isGENIE(RunID):
    # Checks if mcGenerator is GENIE or NEUT
    if int(str(RunID)[1]) == 1:
        return True
    else:
        return False


def getBunchIdx(bunchTimesMean, bunchTimesSigm, pidTime, nSigCutoff=5):
    # We set the time range for each bunch to be:
    # bunchMean +/- nSigCutoff*bunchSig

    # underflow bin
    if 0 < pidTime < T_BUNCH0_CUTOFF:
        return 0
    # Find the corresponding bunch
    for idx, mean in enumerate(bunchTimesMean):
        if (mean - nSigCutoff * bunchTimesSigm[idx] < pidTime < mean + nSigCutoff * bunchTimesSigm[idx]):
            return idx + 1

    # if none found return -1
    return -1


def mcShiftBunchTiming(pidTime, runNum):
    # Bunch timing structure used in GENIE MC prod5C
    bunchTimesMean = [2742., 3324., 3905., 4487., 5069., 5650.]
    bunchTimesSigm = [20., 20., 19.77, 20., 19.28, 20.]
    if runNum > 1:
        bunchTimesMean.extend([6234., 6816.])
        bunchTimesSigm.extend([20., 19.18])

    return getBunchIdx(bunchTimesMean, bunchTimesSigm, pidTime)


def beforeShiftBunchTiming(pidTime, runNum):
    # Bunch timing structure before shift that occurred
    # in early 2011.  THESE ARE SKETCHHHH t.y. 9/5/12
    bunchTimesMean = [2843., 3427., 4008., 4594., 5174., 5757.]
    bunchTimesSigm = [20., 20., 19.77, 20., 19.28, 20.]
    if runNum > 1:
        bunchTimesMean.extend([6334., 6916.])
        bunchTimesSigm.extend([20., 19.18])

    return getBunchIdx(bunchTimesMean, bunchTimesSigm, pidTime)


def afterShiftBunchTiming(pidTime, runNum):
    # Bunch timing structure after shift that occurred
    # in early 2011
    bunchTimesMean = [3017., 3596., 4180., 4764.,
                      5345., 5924., 6506., 7092.]
    bunchTimesSigm = [18.11, 20.51, 20.05, 18.64,
                      19.05, 20.71, 20.18, 19.23]

    return getBunchIdx(bunchTimesMean, bunchTimesSigm, pidTime)


def findBunch(isMC, evtTime, runNum, pidTime):
    # shiftDateUTC corresponds to UTS time when bunch timing shift occurred
    shiftDateUTC = 1293840000
    if not isMC:
        if evtTime < shiftDateUTC:
            return beforeShiftBunchTiming(pidTime, runNum)
        else:
            return afterShiftBunchTiming(pidTime, runNum)
    else:
        # What MC variable gives shift timing status?
        return mcShiftBunchTiming(pidTime, runNum)


def getMaxMomIndices(PIDs, NPIDs, isMC, evtTime, runNum):
    """ This function returns the indices corresponding to the maximum momentum
    tracks in the list of PIDs.  For example, if NPIDs = 3, and PIDs[2] is
    the highest momentum track (sorted bunch by bunch) in bunch 1, then
    this function returns [-1, 2, -1, -1, -1, ..., -1]

    Maximum of 8 bunches per spill + 1 overflow, we need to find the highest
    momentum track in each BUNCH, not in the overall spill.
    """
    maxMom = [-9999999.] * 9
    maxMomIdx = [-1] * 9

    for i in xrange(NPIDs):
        thisPID = PIDs[i]
        reconPos = thisPID.FrontPosition

        # Status of global fit must be good
        if thisPID.Status == 0:
            continue
        # Global fit must have both P0D and Tracker component
        if ((thisPID.NP0Ds < 1) or (thisPID.NTPCs < 1)):
            continue
        # Global vertex must be in P0Dfiducial
        if not inP0DFiducial(reconPos):
            continue

        # If status of fit and both P0D/tracker are used, then reconstructed
        # front momentum is more trustworthy.  So we find the highest
        # momentum track out of all well-fitted tracks.
        reconMom = thisPID.FrontMomentum
        bunchNum = findBunch(isMC, evtTime, runNum, reconPos.T())
        if bunchNum == -1:
            continue
        if reconMom > maxMom[bunchNum]:
            maxMom[bunchNum] = reconMom
            maxMomIdx[bunchNum] = i

    return maxMomIdx


def getPIDsPerBunch(PIDs, NPIDs, isMC, evtTime, runNum):
    """ Returns array[i] with entries corresponding to the number of globalPIDs
    in bunch i.

    Maximum of 8 bunches per spill + 1 overflow and 1 underflow
    """
    pidsPerBunch = [0] * 9
    pidIdxsPerBunch = {}
    for i in range(9):
        pidIdxsPerBunch.setdefault(i, [])

    for i in xrange(NPIDs):
        thisPID = PIDs[i]
        reconPos = thisPID.FrontPosition

        bunchNum = findBunch(isMC, evtTime, runNum, reconPos.T())
        if bunchNum == -1:
            continue

        pidsPerBunch[bunchNum] += 1
        pidIdxsPerBunch[bunchNum].append(i)

    return pidsPerBunch, pidIdxsPerBunch

def sortToBunches(PIDs, NPIDs, isMC, evtTime, runNum):
    """ Combines the functionality of getMaxMomIndices and getPIDsPerBunch.
    Reduces the number of times we loop through the list of PIDs.

    returns three objects of length 9 corresponding the 8 bunches + 1 ovrflw
    maxMomIdx: contains the index of the maximum momentum pid in each bunch
    maxMomNegIdx: contains the index of the maximum momentum neg pid in each bunch
    maxMomPosIdx: contains the index of the maximum momentum pos pid in each bunch
    pidsPerBunch: contains the number of pids in each bunch
    pidIdxsPerBunch: dict containing bunch numbers as keys to list of
    indices of pids in the bunch
    p0dPidsPerBunch: contains the number of pids in each bunch that have 
    nP0D > 0
    p0dPidIdxsPerBunch: dict containing bunch numbers as keys to list of
    indices of p0d pids in the bunch
    """
    maxMom = [-9999999.] * 9
    maxMomNeg = [-9999999.] * 9
    maxMomPos = [-9999999.] * 9
    maxMomIdx = [-1] * 9
    maxMomNegIdx = [-1] * 9
    maxMomPosIdx = [-1] * 9
    pidsPerBunch = [0] * 9
    pidIdxsPerBunch = {}
    p0dPidsPerBunch = [0] * 9
    p0dPidIdxsPerBunch = {}
    for i in range(9):
        pidIdxsPerBunch.setdefault(i, [])
        p0dPidIdxsPerBunch.setdefault(i, [])

    for i in xrange(NPIDs):
        thisPID = PIDs[i]
        reconPos = thisPID.FrontPosition

        bunchNum = findBunch(isMC, evtTime, runNum, reconPos.T())
        if bunchNum == -1:
            continue

        # Status of global fit must be good
        if thisPID.Status == 0:
            continue
        pidsPerBunch[bunchNum] += 1
        pidIdxsPerBunch[bunchNum].append(i)

        # Check if this PID has a P0D component, if not we skip it
        if thisPID.NP0Ds < 1:
            continue
        p0dPidsPerBunch[bunchNum] += 1
        p0dPidIdxsPerBunch[bunchNum].append(i)

        # Check to see if this PID used TPC recon
        if thisPID.NTPCs < 1:
            continue

        # Global vertex of max mom candidate must be in P0Dfiducial
        if not inP0DFiducial(reconPos):
            continue

        # If status of fit and both P0D/tracker are used, then reconstructed
        # front momentum is more trustworthy.  So we find the highest
        # momentum track out of all well-fitted tracks.
        reconMom = thisPID.FrontMomentum
        if reconMom > maxMom[bunchNum]:
            maxMom[bunchNum] = reconMom
            maxMomIdx[bunchNum] = i

        # also sort into HMN, HMP per bunch
        charge = thisPID.Charge
        if charge < 0 and reconMom > maxMomNeg[bunchNum]:
            maxMomNeg[bunchNum] = reconMom
            maxMomNegIdx[bunchNum] = i
        elif charge > 0 and reconMom > maxMomPos[bunchNum]:
            maxMomPos[bunchNum] = reconMom
            maxMomPosIdx[bunchNum] = i

    return (maxMomIdx,
            maxMomNegIdx,
            maxMomPosIdx,
            pidsPerBunch,
            pidIdxsPerBunch,
            p0dPidsPerBunch,
            p0dPidIdxsPerBunch)

def getSandMuIndices(PIDs, NPIDs, isMC, evtTime, runNum):
    # Think the best cut for sand mu is to restrict their front position to the first
    # few layers of the P0D.  Additional requirement that they enter the
    # tracker gets us momentum info.  Looking for 1 track events is probably
    # incorrect after talking to Justyna (T.Y. 6/20/13)

    sandMuIndices = []
    for i in xrange(NPIDs):
        thisPID = PIDs[i]
        reconPos = thisPID.FrontPosition

        # event must fall within bunch time ranges
        bunchNum = findBunch(isMC, evtTime, runNum, reconPos.T())
        if bunchNum <= 0:
            continue
        # Status of global fit must be good
        if thisPID.Status == 0:
            continue
        # Need to use hit info to cut on z-position
        if not inSandRegion(thisPID):
            continue
        # Global fit must have at least a P0D component
        if thisPID.NP0Ds < 1:
            continue
        # Global vertex must be in P0Dsandregion

        sandMuIndices.append(i)

    return sandMuIndices

def isDataFile(file):
    # For water subtraction study we want to get PoT values for:
    # water_data, water_mc, air_data, air_mc
    return 'rdp' in file


def isAirFile(file):
    return 'air' in file or 'RUN3' in file

def getPOTs(filePaths, isToyMC=False):
    # Use this for the water subtraction study
    mcPOTstr = 'POT'
    if isToyMC:
        dataPOTstr = 'ToyMCPOT'
    else:
        dataPOTstr = mcPOTstr

    pot_water_data = 0
    pot_air_data = 0
    pot_water_mc = 0
    pot_air_mc = 0
    for file in filePaths:
        inFile = ROOT.TFile(file)
        t2 = inFile.Get('t2')
        t2.GetEvent(0)

        if isDataFile(file) and not isAirFile(file):
            print file, 'POT is saved as pot_water_data'
            pot_water_data += t2.GetLeaf(dataPOTstr).GetValue()
        elif isDataFile(file) and isAirFile(file):
            print file, 'POT is saved as pot_air_data'
            pot_air_data += t2.GetLeaf(dataPOTstr).GetValue()
        elif not isDataFile(file) and not isAirFile(file):
            print file, 'POT is saved as pot_water_mc'
            pot_water_mc += t2.GetLeaf(mcPOTstr).GetValue()
            if isToyMC:
                pot_water_data += t2.GetLeaf(dataPOTstr).GetValue()
        elif not isDataFile(file) and isAirFile(file):
            print file, 'POT is saved as pot_air_mc'
            pot_air_mc += t2.GetLeaf(mcPOTstr).GetValue()
            if isToyMC:
                pot_air_data += t2.GetLeaf(dataPOTstr).GetValue()
        else:
            sys.exit('WTF is this file?')

        inFile.Close()

    # If we're dealing with a ToyMC sample, we must subtract the MCPOT
    # by the "data" POT since we're treating the ToyMC as separate from
    # the rest of the MC.  The POT counts above are the full POT count
    # for MC.
    if isToyMC:
        pot_water_mc -= pot_water_data
        pot_air_mc -= pot_air_data

    return pot_water_data, pot_air_data, pot_water_mc, pot_air_mc

def getPOTsFromChain(ch):
    # Use this to get POTs from ch containing chained t2 trees
    pot = 0
    # keep track of sand_mc_pot separately
    sand_mc_pot = 0
    for i in range(ch.GetEntries()):
        ch.GetEvent(i)
        fileName = ch.GetFile().GetName()
        thispot = ch.GetLeaf('POT').GetValue()
        # sand mc doesn't keep track of POT correctly in ROOTrackerVtx, must scale
        # by 2.5E17 POT/file
        if 'sand' in fileName:
            if '5F' in fileName:
                sand_mc_pot += thispot * 2.5E17
            else:
                # otherwise assume prod5D neut sand MC
                sand_mc_pot += thispot * 5.0E17
        # POT is set to unity for cosmics
        elif 'Cosmic' in fileName:
            return 1
        else:
            pot += thispot

        print 'TFile:', fileName, ', POT:', pot, ', SandMC POT:', sand_mc_pot

    return pot, sand_mc_pot

##########################################################################
# Binning utils
#
##########################################################################

def scaleBins(hist, binSize):
    """ Scales all bins' bin contents by the appropriate factor for a
    constant bin width.
    Necessary when using variable bin sizes.

    binSize in MeV
    Only use when binning enu!!!!!!!!
    """
    nBins = hist.GetNbinsX()
    for i in range(1, nBins + 1):
        hist.SetBinContent(
            i,
            hist.GetBinContent(i) * binSize / (1000 * hist.GetBinWidth(i)))
        hist.SetBinError(
            i,
            hist.GetBinError(i) * binSize / (1000 * hist.GetBinWidth(i)))
    return

def getEnuBinLowX(meV_bin=50):
    # Returns Energy bins used by beam group
    binLowX = []

    for x in range(0, 200, 4 * meV_bin):
        binLowX.append(float(x))
    for x in range(200, 1200, meV_bin):
        binLowX.append(float(x))
    for x in range(1200, 1500, 2 * meV_bin):
        binLowX.append(float(x))
    for x in range(1500, 3500, 10 * meV_bin):
        binLowX.append(float(x))
    for x in range(3500, 11000, 20 * meV_bin):
        binLowX.append(float(x))

    return binLowX

def getZBinLowX(sand=False, water_bags=False, get_water_bag_z=False):
    """ z-binning dependent on p0dule position.  Taken from raj's code at
    /Raid/users/rajd/T2K/Analysis/NEWPLOTS.C

    water_bags bool determines whether we want to separate the water bags
    from the rest of the water p0dule.  From P0D NIM water bags are 28 mm

    get_water_bag_z: Returns a list of the z-positions at the upstream 
    point of each water bag instead of water_bags and p0dules.
    """
    bin_lowx = []

    # USECAL and upstream of USECAL
    if sand:
        # Most sand muons will be upstream in the P0D so we extend the base
        # N_UPSTREAM_P0DULES for a few additional layers
        bin_lowx.extend([WT_UPSTREAM-Z_EC_WIDTH*i
                         for i in reversed(range(N_UPSTREAM_P0DULES_sand+1))])
    else:
        bin_lowx.extend([WT_UPSTREAM-Z_EC_WIDTH*i
                         for i in reversed(range(N_UPSTREAM_P0DULES+1))])

    # USWT and CWT
    water_bag_z = []
    if water_bags or get_water_bag_z:
        for i in range(N_BGS_USWT):
            water_bag_z.append(bin_lowx[-1]+Z_PD_WIDTH)
            bin_lowx.append(bin_lowx[-1]+Z_PD_WIDTH)
            bin_lowx.append(bin_lowx[-1]+Z_BG_WIDTH)
        # water cover
        bin_lowx.append(bin_lowx[-1]+Z_WC_WIDTH)
        for i in range(N_BGS_CWT):
            water_bag_z.append(bin_lowx[-1]+Z_PD_WIDTH)
            bin_lowx.append(bin_lowx[-1]+Z_PD_WIDTH)
            bin_lowx.append(bin_lowx[-1]+Z_BG_WIDTH)
    else:
        bin_lowx.extend([bin_lowx[-1]+Z_WT_WIDTH*(i+1)
                         for i in range(N_BGS_USWT+N_BGS_CWT)])
    # CWT has one additional P0Dule at downstream end
    bin_lowx.append(bin_lowx[-1]+Z_PD_WIDTH)

    # CECAL
    bin_lowx.extend([bin_lowx[-1]+Z_EC_WIDTH*(i+1)
                     for i in range(N_DOWNSTREAM_P0DULES)])

    return water_bag_z if get_water_bag_z else bin_lowx 

def getBagStartZ():
    """ Returns a list of the z-positions at the upstream point of each
    water bag.  Wrapper around getZBinLowX
    """
    # first get the z-pos of bags and p0dules
    water_bag_z = getZBinLowX(get_water_bag_z=True)
    return water_bag_z

def getXLayers():
    """ Returns a list of tuples containing the beginning and end of the
    X-layers in the P0D
    """
    water_bag_z = getBagStartZ()
    
    return [(bag_z+Z_BG_WIDTH, bag_z+Z_BG_WIDTH+Z_PD_WIDTH/2)
            for bag_z in water_bag_z]
    
def getYLayers():
    """ Returns a list of tuples containing the beginning and end of the
    Y-layers in the P0D
    """
    water_bag_z = getBagStartZ()

    return [(bag_z+Z_BG_WIDTH+Z_PD_WIDTH/2, bag_z+Z_BG_WIDTH+Z_PD_WIDTH)
            for bag_z in water_bag_z]

def getBagXLayers():
    """ Returns a list of tuples containing the beginning and end of both the
    Water-bags and X-layers in the P0D
    """
    water_bag_z = getBagStartZ()

    # each layer is a 2 component tuple with the start of the bag and
    # the end of the water+Xlayer_p0dule
    return [(bag_z, bag_z+Z_BG_WIDTH+Z_PD_WIDTH/2) for bag_z in water_bag_z]

def getLayerCutStr(cut_var, layer):
    """ Returns the ROOT cut string for any z-pos layer in the P0D using
    'cut_var'.  It should be a string of the form '((_0_beg < cut_var &&
    cut_var < _0_end) || (_1_beg < cut_var && cut_var < _1_end) ...)' where
    the _n indicates the layer number.

    cut_var: string that corresponds to the variable on which the cut is
    applied.  Usually this is some z-pos variable.

    layer: string that indicates which p0d-z layer to cut on.  Needs to be
    a key in the layers_dict dictionary, which currently is implemented
    for 'bagx', 'x', and 'y'.
    """
    layers_dict = {'bagx':getBagXLayers(),
                   'x':getXLayers(),
                   'y':getYLayers()}

    layers = layers_dict[layer]
    single_layer_template = '({0[0]} < {1} && {1} < {0[1]})'    
    cut_str = '('
    for a_layer in layers:
        cut_str += single_layer_template.format(a_layer, cut_var)
        cut_str += '||' 

    # strip the last '||' and return
    return cut_str.strip('|&')+')'

def getDirZBinLowX(step=1):
    # Returns CosTheta bins with variable binning
    binLowX = []

    for x in range(-1000, 0, step * 500):
        binLowX.append(x / 1000.)
    for x in range(0, 400, step * 200):
        binLowX.append(x / 1000.)
    for x in range(400, 700, 100 * step):
        binLowX.append(x / 1000.)
    for x in range(700, 900, 50 * step):
        binLowX.append(x / 1000.)
    for x in range(900, 1010, 10 * step):
        binLowX.append(x / 1000.)

    return binLowX

#####################
# File/Job Processing
#####################

def replaceInFile(filePath, modFilePath, pattern, subst):
    try:
        newFile = open(modFilePath, 'w')
        oldFile = open(filePath)
        for line in oldFile:
            newFile.write(line.replace(pattern, subst))

    finally:
        newFile.close
        oldFile.close

def CheckFluxWeightFile(fileList):
    """ We need to modify the flux weighting files depending on which Runs the fileList
    contains.
    """

    # Data files are not flux reweighed.  Don't need to modify parameter file
    if 'rdp' in fileList:
        return

    p0dnumuccanalysis_root = os.environ['P0DNUMUCCANALYSISROOT']
    p0dnumuccanalysis_par = p0dnumuccanalysis_root + '/parameters'
    templateFilePath = p0dnumuccanalysis_par + \
        '/p0dNumuCCAnalysis.template.dat'
    parsFilePath = templateFilePath.replace('template', 'parameters')

    searchStr = '@RUNSPECIFIC@'
    commonStr = '/nd5_tuned11bv3.2_11anom_'
    if 'RUN1' in fileList:
        replaceStr = 'run1' + commonStr + 'run1'
        replaceInFile(templateFilePath, parsFilePath, searchStr, replaceStr)
    elif 'RUN2' in fileList:
        replaceStr = 'run2' + commonStr + 'run2'
        replaceInFile(templateFilePath, parsFilePath, searchStr, replaceStr)
    elif 'RUN3' in fileList:
        replaceStr = 'run3c' + commonStr + 'run3c'
        replaceInFile(templateFilePath, parsFilePath, searchStr, replaceStr)
    else:
        sys.exit('P0DCutUtils.py: Incorrect run in filename ' + fileList)

def CheckSubmittedJob(jobFilePath):
    # Quick workaround to ensure taht the correct parameters file is being read in
    # for each cluster submitted process.  We need to make sure that the RunAnalysis.exe
    # has read in the parameters file before attempting to modify it for the next
    # submitted job.  The parameters file contains the flux ratio filepath, which
    # is different depending on the processed run.
    while True:
        if os.path.isfile(jobFilePath):
            if 'Running on file' in open(jobFilePath).read():
                break

        time.sleep(30)

#####################
# Prompts and Options
#####################
# useful prompts for plotting options

def PromptOutputDir(query):
    dir = raw_input(query)
    if not os.path.isdir(dir):
        print 'P0DCutUtils.py: No such outfile dir! ' + dir

        dir = PromptOutputDir(query)

    return dir + '/'


def PromptYesNo(query):
    boolStr = raw_input(query)
    if boolStr.capitalize()[0] == 'Y':
        bool = True
    elif boolStr.capitalize()[0] == 'N':
        bool = False
    else:
        bool = PromptYesNo(query)

    return bool


def PromptVariableToPlot():
    varStr = raw_input('What variable should I plot? ')
    return SearchVar(varStr.upper())


def SearchVar(genericVar):
    # THis function returns the correct variable to be plotted for the highLevelAnalysis
    # tools and for my old analyses (WaterSubtraction) since they are different
    # Returns a tuple of the form genericVar, ('highLevelVar', 'myAnalVar',
    # 'units')
    variableDict = {'FMOMENTUM': ('selmu_mom/1e3',
                                  'GlobalFMom/1e3',
                                  '[GeV]'),
                    'COSTHETA': ('selmu_dir[2]',
                                 'GlobalFDir.CosTheta()',
                                 ''),
                    'CHARGE': ('selmu_charge',
                               'Charge',
                               '')}

    if genericVar in variableDict:
        return genericVar, variableDict[genericVar]
    else:
        print 'P0DCutUtils.SearchVar(): Undefined variable! ' + genericVar
        return PromptVariableToPlot()


def PromptRangeForVar():
    try:
        nBins = int(raw_input('nBins? '))
        while nBins <= 0:
            nBins = int(raw_input('nBins? '))
        binLow = float(raw_input('binLow? '))
        binHigh = float(raw_input('binHigh? '))
        while binLow >= binHigh:
            binLow = float(raw_input('binLow? '))
            binHigh = float(raw_input('binHigh? '))

        binWidth = (binHigh - binLow)/nBins
        return array('d', [binLow+i*binWidth for i in range(nBins+1)])
    except:
        PromptRangeForVar()


def PromptSaveOption():
    usrin = raw_input(
        'How should output plots be saved to ROOT file? (UPDATE or RECREATE) ')
    if usrin.capitalize()[0] == 'R':
        return 'RECREATE'
    else:
        return 'UPDATE'

def prompt(query, check_list):
    # Prompt and check that what the user entered is in a list, check_list
    usr_in = raw_input(query)

    if usr_in == '':
        print 'Setting to default value, the first item in list:', check_list[0]
        return check_list[0]
    elif usr_in == '!':
        print 'Canceling... nothing set'
        return

    check = usr_in in check_list
    if not check:
        return prompt(query, check_list)

    return usr_in

################################## TPC-P0D matching stuff ################


def projectIntoP0D(tpc1_front_pos, tpc1_front_dir):
    """ Projects position of tpc1-track into P0D using the tpc1_front_dir

    returns: TVector of projected position
    """

    projected_z_plane = -1100.

    delta_z = tpc1_front_pos.Z() - projected_z_plane

    # Calculate delta-y, delta-x assuming linear projection
    delta_y = tpc1_front_dir.Y() / tpc1_front_dir.Z() * delta_z
    delta_x = tpc1_front_dir.X() / tpc1_front_dir.Z() * delta_z

    projected_y = tpc1_front_pos.Y() - delta_y
    projected_x = tpc1_front_pos.X() - delta_x

    return ROOT.TVector3(projected_x, projected_y, projected_z_plane)


def convertToDict(chain, **vars_dict):
    """ Convert a tree into a dict with 'runID/eventID/var_name' as keys.
    This will reduce memory usage by half compared to using a nested dict
    structure when loading dict into memory.

    e.g.: tree contains RunID:EventID:Some_Var.  Then the dict will be of
    the form: {'RunNum/EventNum/var_1_name':var_1_val, ... }

    vars_dict: dicts of the form {var_1_name:var_1_type, ...}.
    Currently only works for simple types.
    """

    # Keep track of the arrays of variables matched to each tbranch
    # var_arrays: {var_name: array(var_type, [0])
    var_arrays = {}
    for a_var_name, a_var_type in vars_dict.iteritems():
        # Create an array of type a_var_type to store the variable
        var_array = array(a_var_type, [0])
        chain.SetBranchAddress(a_var_name, var_array)

        var_arrays[a_var_name] = var_array

    converted = {}
    key_base = '{run}\{event}\{var}'
    # loop through chain and store in dictionary
    for e in range(chain.GetEntries()):
        chain.GetEntry(e)
        run_id = chain.RunID
        event_id = chain.EventID

        # convert the array to a python value, and store in 1D dictionary
        for name, arr in var_arrays.iteritems():
            key_name = key_base.format(run = run_id,
                                       event = event_id,
                                       var = name)
            converted[key_name] = arr[0]

    # print converted
    return converted


def convertToNestedDict(chain, **vars_dict):
    """ Convert a tree into a dict-of-dicts with runID and eventID as keys.

    e.g.: tree contains RunID:EventID:Some_Var.  Then the dict will be of
    the form: {RunNum:{EventNum: {var_1_name:var_1_val, ... }, ... }, ... }

    vars_dict: dicts of the form {var_1_name:var_1_type, ...}.
    Currently only works for simple types.
    """

    # Keep track of the arrays of variables matched to each tbranch
    # var_arrays: {var_name: array(var_type, [0])
    var_arrays = {}
    for a_var_name, a_var_type in vars_dict.iteritems():
        # Create an array of type a_var_type to store the variable
        var_array = array(a_var_type, [0])
        chain.SetBranchAddress(a_var_name, var_array)

        var_arrays[a_var_name] = var_array

    # Autovivification class makes nested dictionaries easier to
    # handle.  We can do blah[1][2][3]=4
    converted = autovivification.AutoVivification()
    # loop through chain and store in dictionary
    for e in range(chain.GetEntries()):
        chain.GetEntry(e)
        run_id = chain.RunID
        event_id = chain.EventID

        # convert the array to a python value, and store in
        # the autovivified dictionary (see google)
        for name, arr in var_arrays.iteritems():
            converted[run_id][event_id][name] = arr[0]

    # print converted
    return converted
