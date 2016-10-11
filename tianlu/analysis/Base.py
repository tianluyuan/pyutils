"""
Base.py: Set all variables we want to extract from oaAnalysis in __init__
This way all derived classes such as GetEvents or GetSandEvents will have
the same initialization.
Also the class includes some useful utility methods that can then be called
within child classes.  The methods that derive from userAnalysisBase can
be overwritten in the child class.
"""

# import array module so can give ROOT its pointers
from array import array
# import the ROOT module
import ROOT
# import the userAnalysisBase class from userAnalysisBase.py
from userAnalysisBase import userAnalysisBase
# import P0D cut utils
import sys
import os
import glob
import marshal
import time
from numpy import math
from tianlu import directories, rootutils
from tianlu.analysis import utils, P0DCutUtils
from exceptions import AttributeError
from lockfile import FileLock


__author__ = 'Tianlu Yuan'
__email__ = 'tianlu.yuan [at] colorado.edu'


class P0DNumuAnalysisBase(userAnalysisBase, object):
    """ Derives from userAnalysisBase, and is called from runAnalyses.py
    class name must be the same as the filename without the .py
    """
    def __init__(self, options):
        """ Specifies further trees to be added to oaAnalysisModuleTChains.
        """
        # Call __init__ in base class
        userAnalysisBase.__init__(self, options)
        # Make TChains of the oaAnalysis trees to be read, and append to i
        # oaAnalysisModuleTChains
        self.BasicHeader = ROOT.TChain('HeaderDir/BasicHeader')
        self.BasicDataQuality = ROOT.TChain('HeaderDir/BasicDataQuality')
        self.BeamSummaryData = ROOT.TChain('HeaderDir/BeamSummaryData')
        self.Global = ROOT.TChain('ReconDir/Global')
        self.P0D = ROOT.TChain('ReconDir/P0D')
        self.Trajectories = ROOT.TChain('TruthDir/Trajectories')
        self.Vertices = ROOT.TChain('TruthDir/Vertices')
        self.GRooTrackerVtx = ROOT.TChain('TruthDir/GRooTrackerVtx')
        self.NRooTrackerVtx = ROOT.TChain('TruthDir/NRooTrackerVtx')

        self.oaAnalysisModuleTChains.append(self.BasicHeader)
        self.oaAnalysisModuleTChains.append(self.BasicDataQuality)
        self.oaAnalysisModuleTChains.append(self.BeamSummaryData)
        self.oaAnalysisModuleTChains.append(self.Global)
        self.oaAnalysisModuleTChains.append(self.P0D)
        self.oaAnalysisModuleTChains.append(self.Trajectories)
        self.oaAnalysisModuleTChains.append(self.Vertices)
        self.oaAnalysisModuleTChains.append(self.GRooTrackerVtx)
        self.oaAnalysisModuleTChains.append(self.NRooTrackerVtx)

        # Counter to keep track MC trees when adding MC PoTs
        self.NextPOTTreeNumber = 0

        # Flux tuning histos from flux_version 13av1.0, nd5 files
        # nd5 is the ND280 Basket, the flux files of which I'm
        # assuming is similar enough to be used for ND280 Magnet files.
        # For the tuning, use "*_fine.root" and choose tuned13a* histos,
        # which are the histos with total tuning applied.  The "real*"
        # histos in the root file only have pbeam tuning applied.
        # See http://www.t2k.org/beam/NuFlux/FluxRelease/11arelease

        # TODO: Fix so that we pick the corresponding flux file for each run
        # not the avg.
        fluxFilePath = utils.flux_ratio_file()
        fluxFile = ROOT.TFile(fluxFilePath)
        # We want the ratio between the flux using real P-beam profile and 11a
        # nominal.
        self.fluxRatio_numu = fluxFile.Get('enu_nd5_tuned13a_numu_ratio')
        self.fluxRatio_numub = fluxFile.Get('enu_nd5_tuned13a_numub_ratio')
        self.fluxRatio_nue = fluxFile.Get('enu_nd5_tuned13a_nue_ratio')
        self.fluxRatio_nueb = fluxFile.Get('enu_nd5_tuned13a_nueb_ratio')

        # Python needs to give up ownership of fluxFile, otherwise it is
        # closed when __init__ is finished
        ROOT.SetOwnership(fluxFile, False)

        # self.listOfHistosAndOtherObjects.extend([self.fluxRatio_numu, self.fluxRatio_numub,
        # self.fluxRatio_nue, self.fluxRatio_nueb])

        # Dictionary to keep track of the MCM-Time-Since-Busy as derived from
        # reco
        self._mcm_time_dict = None
        # size of arrays to store p0d-hit-info
        self._big_array_size = 2000
        # list of runs that we need to flag as bad
        # Ignore run 8359 (part of Run3b), seems to be undergoing P0D-LI-calib.
        self._bad_run_list = [8359]
        self.tianluFlag = False


    def BookHistogramsAndAnalysisSteps(self):
        """ Derived method from userAnalysisBase
        Bookkeeping
        """

        # Initiates an eCorr object
        self.eCorr = ROOT.P0DELossCorr()

        # Define variables for output tree. As pointer don't exist in python we
        # have to declare an array and then pass this as the branch address as
        # ROOT expects a pointer.
        #
        # Truth
        self._FPos = ROOT.TLorentzVector()
        self._FMom = ROOT.TLorentzVector()
        self._VtxPos = ROOT.TLorentzVector()
        self._TrueVtxInP0DFiducial = array("i", [0])
        self._TruePDG = array("i", [0])
        self._NuMom = ROOT.TLorentzVector()
        self._NuPDG = array("i", [0])
        self._TgtMom = ROOT.TLorentzVector()
        self._TargetPDG = array("i", [0])
        self._Reaction = ROOT.TString()
        self._ReactionInt = array("i", [0])
        self._Eff = array("d", [0.0])
        self._Pur = array("d", [0.0])
        # Global
        self._GlobalFPos = ROOT.TLorentzVector()
        self._GlobalFDir = ROOT.TVector3()
        self._GlobalFMom = array("d", [0.0])
        self._GlobalCone = ROOT.TVector3()
        self._Charge = array("d", [0.0])
        self._BunchNum = array("i", [0])
        self._PIDsInBunch = array("i", [0])
        self._P0DPIDsInBunch = array("i", [0])
        self._nTPCs = array("i", [0])
        self._nP0Ds = array("i", [0])
        self._nPDGs = array("i", [0])
        self._PDGs = array("i", [0] * 10)
        self._PDGWeights = array("d", [0.0] * 10)
        self._isForward = array("i", [0])
        self._goodGlobalTrack = array("i", [0])
        # TPCObject from GlobalPID info
        self._TPC1Used = array("i", [0])
        self._TPC1NNodes = array("i", [0])
        self._TPC1Len = array("d", [0.0])
        self._TPC1FPos = ROOT.TLorentzVector()
        self._TPC1FMom = array("d", [0.0])
        self._TPC1FDir = ROOT.TVector3()
        self._TPC1Charge = array("d", [0.0])
        self._TPC1Chi2 = array("d", [0.0])
        # P0DObject from GlobalPID info
        self._P0DLen = array("d", [0.0])
        self._P0DPosDiff = array("d", [0.0])
        self._P0DELoss = array("d", [0.0])
        self._P0DRevELoss = array("d", [0.0])
        self._P0DFPos = ROOT.TLorentzVector()
        self._P0DBPos = ROOT.TLorentzVector()
        self._P0DFMom = array("d", [0.0])
        self._P0DBMom = array("d", [0.0])
        self._P0DCone = ROOT.TVector3()
        self._P0DFDir = ROOT.TVector3()
        # Truth-Global Matching
        # Integer that returns true if track found using GlobalPID matches
        # true muon track.
        self._trueGlobalMatched = array("i", [0])
        self._vtxGlobalMatched = array("i", [0])
        # TPC-P0D Matching
        self._P0DGoing = array("i", [0])
        self._P0DActive = array("i", [0])
        self._nGoodHits = array("i", [0])
        self._MCMTimeSinceBusy = array("i", [0])
        self._nP0DHits = array("i", [0])
        self._P0DReconTrackMatch = array("i", [0])
        self._P0DReconShowerMatch = array("i", [0])
        self._P0DGlobalTrackMatch = array("i", [0])
        self._P0DGlobalShowerMatch = array("i", [0])
        self._P0DHitChanID = array("I", [0] * self._big_array_size)
        self._P0DHitTimes = array("f", [0.0] * self._big_array_size)
        self._P0DAdjHitTimes = array("f", [0.0] * self._big_array_size)
        self._P0DAdjHitTimes2 = array("f", [0.0] * self._big_array_size)
        self._P0DHitCharges = array("f", [0.0] * self._big_array_size)
        self._TPC1AdjTime = array("f", [0])
        # P0DRecon Info, from matched GlobalPID p0d object
        self._P0DRecoNMichels = array("i", [0])
        self._P0DRecoTrackFound = array("i", [0])
        self._P0DRecoShowerFound = array("i", [0])
        self._P0DRecoELoss = array("d", [0.0])
        self._P0DRecoRevELoss = array("d", [0.0])
        # Header info
        self._IsMC = array("i", [0])
        self._IsSandMC = array("i", [0])
        self._ToyMC = array("i", [0])
        self._P0DWaterStatus = array("i", [0])
        self._RunID = array("i", [0])
        self._EventID = array("i", [0])
        self._SubrunID = array("i", [0])
        self._EventTime = array("l", [0])
        self._FGDCosmic = array("i", [0])
        self._TripTCosmic = array("i", [0])
        self._P0DFlag = array("i", [0])
        self._P0DECALFlag = array("i", [0])
        self._FGDFlag = array("i", [0])
        self._TPCFlag = array("i", [0])
        self._ECALFlag = array("i", [0])
        self._DSECALFlag = array("i", [0])
        # Flux tuning weight
        self._FluxRatio = array("d", [0.0])

        # Create output tree and set branch addresses.
        self.outputTree = ROOT.TTree("t1", "Yiggdrasil")
        # Header Branches
        self.outputTree.Branch("RunID", self._RunID, "RunID/I")
        self.outputTree.Branch("EventID", self._EventID, "EventID/I")
        self.outputTree.Branch("SubrunID", self._SubrunID, "SubrunID/I")
        self.outputTree.Branch("EventTime", self._EventTime, "EventTime/I")
        self.outputTree.Branch("IsMC", self._IsMC, "IsMC/B")
        self.outputTree.Branch("IsSandMC", self._IsSandMC, "IsSandMC/B")
        self.outputTree.Branch("ToyMC", self._ToyMC, "ToyMC/I")
        self.outputTree.Branch(
            "P0DWaterStatus",
            self._P0DWaterStatus,
            "P0DWaterStatus/B")
        self.outputTree.Branch("FGDCosmic", self._FGDCosmic, "FGDCosmic/I")
        self.outputTree.Branch(
            "TripTCosmic",
            self._TripTCosmic,
            "TripTCosmic/I")
        self.outputTree.Branch("P0DFlag", self._P0DFlag, "P0DFlag/I")
        self.outputTree.Branch(
            "P0DECALFlag",
            self._P0DECALFlag,
            "P0DECALFlag/I")
        self.outputTree.Branch("FGDFlag", self._FGDFlag, "FGDFlag/I")
        self.outputTree.Branch("TPCFlag", self._TPCFlag, "TPCFlag/I")
        self.outputTree.Branch("ECALFlag", self._ECALFlag, "ECALFlag/I")
        self.outputTree.Branch("DSECALFlag", self._DSECALFlag, "DSECALFlag/I")
        # Truth Branches
        self.outputTree.Branch(
            "InitPos",
            "TLorentzVector",
            ROOT.AddressOf(self._FPos))
        self.outputTree.Branch(
            "InitMom",
            "TLorentzVector",
            ROOT.AddressOf(self._FMom))
        self.outputTree.Branch(
            "VtxPos",
            "TLorentzVector",
            ROOT.AddressOf(self._VtxPos))
        self.outputTree.Branch(
            "TrueVtxInP0DFiducial",
            self._TrueVtxInP0DFiducial,
            "TrueVtxInP0DFiducial/B")
        self.outputTree.Branch("TruePDG", self._TruePDG, "TruePDG/I")
        self.outputTree.Branch(
            "NeutrinoMomentum",
            "TLorentzVector",
            ROOT.AddressOf(self._NuMom))
        self.outputTree.Branch("NeutrinoPDG", self._NuPDG, "NeutrinoPDG/I")
        self.outputTree.Branch(
            "TargetMomentum",
            "TLorentzVector",
            ROOT.AddressOf(self._TgtMom))
        self.outputTree.Branch("TargetPDG", self._TargetPDG, "TargetPDG/I")
        self.outputTree.Branch(
            "Reaction",
            "TString",
            ROOT.AddressOf(self._Reaction))
        self.outputTree.Branch(
            "ReactionInt",
            self._ReactionInt,
            "ReactionInt/I")
        self.outputTree.Branch("Eff", self._Eff, "Eff/D")
        self.outputTree.Branch("Pur", self._Pur, "Pur/D")
        # Global Branches
        self.outputTree.Branch("Charge", self._Charge, "GlobalCharge/D")
        self.outputTree.Branch(
            "GlobalFPos",
            "TLorentzVector",
            ROOT.AddressOf(self._GlobalFPos))
        self.outputTree.Branch(
            "GlobalFDir",
            "TVector3",
            ROOT.AddressOf(self._GlobalFDir))
        self.outputTree.Branch(
            "GlobalFMom",
            self._GlobalFMom,
            "GlobalFrontMomentum/D")
        self.outputTree.Branch(
            "GlobalCone",
            "TVector3",
            ROOT.AddressOf(self._GlobalCone))
        self.outputTree.Branch("BunchNum", self._BunchNum, "BunchNum/I")
        self.outputTree.Branch(
            "PIDsInBunch",
            self._PIDsInBunch,
            "PIDsInBunch/I")
        self.outputTree.Branch(
            "P0DPIDsInBunch",
            self._P0DPIDsInBunch,
            "P0DPIDsInBunch/I")
        self.outputTree.Branch("nTPCs", self._nTPCs, "nTPCs/I")
        self.outputTree.Branch("nP0Ds", self._nP0Ds, "nP0Ds/I")
        self.outputTree.Branch("nPDGs", self._nPDGs, "nPDGs/I")
        self.outputTree.Branch("PDGCodes", self._PDGs, "PDGs[nPDGs]/I")
        self.outputTree.Branch(
            "PDGWeights",
            self._PDGWeights,
            "PDGWeights[nPDGs]/D")
        self.outputTree.Branch("isForward", self._isForward, "isForward/I")
        self.outputTree.Branch(
            "goodGlobalTrack",
            self._goodGlobalTrack,
            "goodGlobalTrack/I")
        # TPC1 from GlobalPID
        self.outputTree.Branch("TPC1Used", self._TPC1Used, "TPC1Used/I")
        self.outputTree.Branch("TPC1NNodes", self._TPC1NNodes, "TPC1NNodes/I")
        self.outputTree.Branch("TPC1Len", self._TPC1Len, "TPC1Length/D")
        self.outputTree.Branch(
            "TPC1FPos",
            "TLorentzVector",
            ROOT.AddressOf(self._TPC1FPos))
        self.outputTree.Branch(
            "TPC1FMom",
            self._TPC1FMom,
            "TPC1FrontMomentum/D")
        self.outputTree.Branch(
            "TPC1FDir",
            "TVector3",
            ROOT.AddressOf(self._TPC1FDir))
        self.outputTree.Branch("TPC1Charge", self._TPC1Charge, "TPC1Charge/D")
        self.outputTree.Branch("TPC1Chi2", self._TPC1Chi2, "TPC1Chi2/D")
        # P0D from GlobalPID
        self.outputTree.Branch("P0DLen", self._P0DLen, "P0DLength/D")
        self.outputTree.Branch(
            "P0DPosDiff",
            self._P0DPosDiff,
            "P0DPositionDiff/D")
        self.outputTree.Branch("P0DELoss", self._P0DELoss, "P0DELoss/D")
        self.outputTree.Branch(
            "P0DRevELoss",
            self._P0DRevELoss,
            "P0DRevELoss/D")
        self.outputTree.Branch(
            "P0DFPos",
            "TLorentzVector",
            ROOT.AddressOf(self._P0DFPos))
        self.outputTree.Branch(
            "P0DBPos",
            "TLorentzVector",
            ROOT.AddressOf(self._P0DBPos))
        self.outputTree.Branch("P0DFMom", self._P0DFMom, "P0DFrontMomentum/D")
        self.outputTree.Branch("P0DBMom", self._P0DBMom, "P0DBackMomentum/D")
        self.outputTree.Branch(
            "P0DCone",
            "TVector3",
            ROOT.AddressOf(self._P0DCone))
        self.outputTree.Branch(
            "P0DFDir",
            "TVector3",
            ROOT.AddressOf(self._P0DFDir))
        # Truth-Global Matching
        self.outputTree.Branch(
            "trueGlobalMatch",
            self._trueGlobalMatched,
            "trueGlobalMatch/B")
        self.outputTree.Branch(
            "vtxGlobalMatch",
            self._vtxGlobalMatched,
            "vtxGlobalMatch/B")
        # TPC-P0D Matching
        self.outputTree.Branch("P0DGoing", self._P0DGoing, "P0DGoing/I")
        self.outputTree.Branch("P0DActive", self._P0DActive, "P0DActive/I")
        self.outputTree.Branch("nGoodHits", self._nGoodHits, "nGoodHits/I")
        self.outputTree.Branch(
            "MCMTimeSinceBusy",
            self._MCMTimeSinceBusy,
            "MCMTimeSinceBusy/I")
        self.outputTree.Branch("nP0DHits", self._nP0DHits, "nP0DHits/I")
        self.outputTree.Branch(
            "P0DReconTrackMatch",
            self._P0DReconTrackMatch,
            "P0DReconTrackMatch/B")
        self.outputTree.Branch(
            "P0DReconShowerMatch",
            self._P0DReconShowerMatch,
            "P0DReconShowerMatch/B")
        self.outputTree.Branch(
            "P0DGlobalTrackMatch",
            self._P0DGlobalTrackMatch,
            "P0DGlobalTrackMatch/B")
        self.outputTree.Branch(
            "P0DGlobalShowerMatch",
            self._P0DGlobalShowerMatch,
            "P0DGlobalShowerMatch/B")
        self.outputTree.Branch(
            "P0DHitChanID",
            self._P0DHitChanID,
            "P0DHitChanID[nP0DHits]/s")
        self.outputTree.Branch(
            "P0DHitTimes",
            self._P0DHitTimes,
            "P0DHitTimes[nP0DHits]/F")
        self.outputTree.Branch(
            "P0DAdjHitTimes",
            self._P0DAdjHitTimes,
            "P0DAdjHitTimes[nP0DHits]/F")
        self.outputTree.Branch(
            "P0DAdjHitTimes2",
            self._P0DAdjHitTimes2,
            "P0DAdjHitTimes2[nP0DHits]/F")
        self.outputTree.Branch(
            "P0DHitCharges",
            self._P0DHitCharges,
            "P0DHitCharges[nP0DHits]/F")
        self.outputTree.Branch(
            "TPC1AdjTime",
            self._TPC1AdjTime,
            "TPC1AdjTime/F")
        # P0DReco Info taken from P0DRecon tree as opposed to PID.P0D
        self.outputTree.Branch(
            "P0DRecoNMichels",
            self._P0DRecoNMichels,
            "P0DRecoNMichels/I")
        self.outputTree.Branch(
            "P0DRecoTrackFound",
            self._P0DRecoTrackFound,
            "P0DRecoTrackFound/I")
        self.outputTree.Branch(
            "P0DRecoShowerFound",
            self._P0DRecoShowerFound,
            "P0DRecoShowerFound/I")
        self.outputTree.Branch(
            "P0DRecoELoss",
            self._P0DRecoELoss,
            "P0DRecoELoss/D")
        self.outputTree.Branch(
            "P0DRecoRevELoss",
            self._P0DRecoRevELoss,
            "P0DRecoRevELoss/D")
        # Flux tuning weight
        self.outputTree.Branch("FluxRatio", self._FluxRatio, "FluxRatio/D")

        ###
        # Define variables for counterTree to store POT, num passing cuts, etc
        # info...
        #
        self._POT = array("d", [0.0])
        self._ToyMCPOT = array("d", [0.0])

        # Create Header info tree
        self.counterTree = ROOT.TTree("t2", "Bonsai")

        # Header Branch Consisting of
        self.counterTree.Branch("POT", self._POT, "POT/D")
        self.counterTree.Branch("ToyMCPOT", self._ToyMCPOT, "ToyMCPOT/D")

        #########
        # Define variables for potTree to histogram event rates etc.
        #
        # In order to histogram event rates, we need the PCT5perSpill,
        # eventTimes for each spill, and a flag for the number of good
        # events that pass our cuts in each spill.  Then
        # eventRate = (spillTimes weighted by nEvents)/(spillTimes weighted by
        # PCT5perSpill)
        self._PCT5perSpill = array("d", [0.0])
        self._SpillTime = array("l", [0])
        self._PIDsPerBunchArray = array("i", [0] * 9)
        self._SelectedPIDsPerBunchArray = array("i", [0] * 9)
        self._nBunchesWPIDPerSpill = array("i", [0])
        self._nEventsPerSpill = array("i", [0])
        self._nTPC1UsedPerSpill = array("i", [0])
        self._nSinglePIDBunchPerSpill = array("i", [0])
        self._nSinglePIDwnoEcalPerSpill = array("i", [0])
        self._nSinglePIDwnoP0DEcalPerSpill = array("i", [0])

        # Create Header info tree
        self.potTree = ROOT.TTree("t3", "Protons")

        # Header Branch Consisting of
        self.potTree.Branch(
            "PCT5perSpill",
            self._PCT5perSpill,
            "PCT5perSpill/D")
        self.potTree.Branch("SpillTime", self._SpillTime, "SpillTime/I")
        self.potTree.Branch(
            "PIDsPerBunchArray",
            self._PIDsPerBunchArray,
            "PIDsPerBunchArray[9]/I")
        self.potTree.Branch(
            "SelectedPIDsPerBunchArray",
            self._SelectedPIDsPerBunchArray,
            "SelectedPIDsPerBunchArray[9]/I")
        self.potTree.Branch(
            "nBunchesWPIDPerSpill",
            self._nBunchesWPIDPerSpill,
            "nBunchesWPIDPerSpill/I")
        self.potTree.Branch(
            "nEventsPerSpill",
            self._nEventsPerSpill,
            "nEventsPerSpill/I")
        self.potTree.Branch(
            "nTPC1UsedPerSpill",
            self._nTPC1UsedPerSpill,
            "nTPC1UsedPerSpill/I")
        self.potTree.Branch(
            "nSinglePIDBunchPerSpill",
            self._nSinglePIDBunchPerSpill,
            "nSinglePIDBunchPerSpill/I")
        self.potTree.Branch(
            "nSinglePIDwnoEcalPerSpill",
            self._nSinglePIDwnoEcalPerSpill,
            "nSinglePIDwnoEcalPerSpill/I")
        self.potTree.Branch(
            "nSinglePIDwnoP0DEcalPerSpill",
            self._nSinglePIDwnoP0DEcalPerSpill,
            "nSinglePIDwnoP0DEcalPerSpill/I")

        ########
        # Truth tree
        self._Generator = ROOT.TString()
        self._InitMom = ROOT.TLorentzVector()
        self._InitPos = ROOT.TLorentzVector()
        self._VertPos = ROOT.TLorentzVector()
        self._NeutrinoMomentum = ROOT.TLorentzVector()
        self._TargetMomentum = ROOT.TLorentzVector()
        self._NeutrinoPDG = array("i", [0])
        self._TargetPDG = array("i", [0])
        self._TrajPDG = array("i", [0])
        self._ID = array("i", [0])
        self._ReactionCode = ROOT.TString()
        self._ReactionConverted = array("i", [0])

        # Create tree
        self.truthTree = ROOT.TTree("t4", "Truth Selection")

        # Create branches
        self.truthTree.Branch(
            "RunID",
            self._RunID,
            "RunID/I")
        self.truthTree.Branch(
            "EventID",
            self._EventID,
            "EventID/I")
        self.truthTree.Branch(
            "SubrunID",
            self._SubrunID,
            "SubrunID/I")
        self.truthTree.Branch(
            "Generator",
            "TString",
            ROOT.AddressOf(self._Generator))
        self.truthTree.Branch(
            "InitPos",
            "TLorentzVector",
            ROOT.AddressOf(self._InitPos))
        self.truthTree.Branch(
            "InitMom",
            "TLorentzVector",
            ROOT.AddressOf(self._InitMom))
        self.truthTree.Branch(
            "VtxPos",
            "TLorentzVector",
            ROOT.AddressOf(self._VertPos))
        self.truthTree.Branch(
            "NeutrinoMomentum",
            "TLorentzVector",
            ROOT.AddressOf(self._NeutrinoMomentum))
        self.truthTree.Branch(
            "TargetMomentum",
            "TLorentzVector",
            ROOT.AddressOf(self._TargetMomentum))
        self.truthTree.Branch(
            "NeutrinoPDG",
            self._NeutrinoPDG,
            "NeutrinoPDG/I")
        self.truthTree.Branch(
            "TargetPDG",
            self._TargetPDG,
            "TargetPDG/I")
        self.truthTree.Branch(
            "TrajPDG",
            self._TrajPDG,
            "TrajPDG/I")
        self.truthTree.Branch(
            "ID",
            self._ID,
            "ID/I")
        self.truthTree.Branch(
            "ReactionCode",
            "TString",
            ROOT.AddressOf(self._ReactionCode))
        self.truthTree.Branch(
            "ReactionConverted",
            self._ReactionConverted,
            "ReactionConverted/I")
 
        ########
        # Define a tree to store variables from GlobalRecon::GlobalVertex class
        # Want to see if vertices with multiple tracks have improved resolution
        # Compare with TrueVtxPos
        #
        self._VertexPos = ROOT.TLorentzVector()
        self._NConstituents = array("i", [0])
        self._TrueVtxPos = ROOT.TLorentzVector()

        # Create Vertex info tree
        self.verticesTree = ROOT.TTree("t5", "Vtxs")

        # Vertex Branch Consisting of
        self.verticesTree.Branch(
            "VertexPos",
            "TLorentzVector",
            ROOT.AddressOf(self._VertexPos))
        self.verticesTree.Branch(
            "NConstituents",
            self._NConstituents,
            "NConstituents/I")
        self.verticesTree.Branch(
            "TrueVtxPos",
            "TLorentzVector",
            ROOT.AddressOf(self._TrueVtxPos))

        #         self.MyHist = ROOT.TH1D('myHist', '', 100, -1500,1500)


    def WriteResultsToOutputFile(self):
        """
        Derived method from userAnalysisBase
        """
        if self.rootOutput:
            self.rootOutput.cd()
            for histo in self.listOfHistosAndOtherObjects:
                histo.Write()
            self.outputTree.Write()
            self.counterTree.Fill()
            self.counterTree.Write()
            self.potTree.Write()
            self.truthTree.Write()
            self.verticesTree.Write()
            summary = ROOT.TObjString(self.GetAnalysisCutResults())
            summary.Write('AnalysisSummary')
            self.rootOutput.Close()
            print "Output saved in " + self.options.output
        else:
            print "Unable to find output file. Tree not saved!"
        # self.MyHist.Draw()
        # raw_input()


    def FillPOTTree(self):
        """ Fill and reset POT counting tree
        """
        self.potTree.Fill()

        # need to reset counters to 0
        # self._PCT5perSpill[0] = 0
        # self._SpillTime[0] = 0
        for i in range(len(self._SelectedPIDsPerBunchArray)):
            self._SelectedPIDsPerBunchArray[i] = 0
        self._nBunchesWPIDPerSpill[0] = 0
        self._nEventsPerSpill[0] = 0
        self._nTPC1UsedPerSpill[0] = 0
        self._nSinglePIDBunchPerSpill[0] = 0
        self._nSinglePIDwnoEcalPerSpill[0] = 0
        self._nSinglePIDwnoP0DEcalPerSpill[0] = 0


    def SetupEventsAnalysis(self):
        """ This is used to setup the variables that will be used to store
        all the relevant info
        """
        self.BeamData = self.BeamSummaryData.BeamSummaryData.At(0)
        self.PIDs = self.Global.PIDs

        # Check MC status quick hack pyroot doesn't handle bools well.
        # self._IsMC[0] = runID > 100000
        self._IsMC[0] = bool(self.BasicHeader.IsMC)
#         print self.BasicHeader.P0DWaterStatus
#         print bool(self.BasicHeader.P0DWaterStatus)
        self._P0DWaterStatus[0] = bool(self.BasicHeader.P0DWaterStatus)
        self._RunID[0] = self.BasicHeader.RunID
        self._SubrunID[0] = self.BasicHeader.SubrunID
        self._EventID[0] = self.BasicHeader.EventID
        self._EventTime[0] = self.BasicHeader.EventTime
        self._FGDCosmic[0] = self.BasicHeader.FGDCosmicEvent
        self._TripTCosmic[0] = self.BasicHeader.TripTCosmicEvent

        # check runID volume integer to see if mcp sand
        self.iscosmic = self._FGDCosmic[0] or self._TripTCosmic[0]
        if self._IsMC[0] and not self.iscosmic:
            self._IsSandMC[0] = bool(str(self._RunID[0])[4] == '7')

        # Detector flags, all should be 0 for those events that pass
        # DataQuality cuts
        self._P0DFlag[0] = self.BasicDataQuality.P0DFlag
        self._P0DECALFlag[0] = self.BasicDataQuality.P0DECALFlag
        self._FGDFlag[0] = self.BasicDataQuality.FGDFlag
        self._TPCFlag[0] = self.BasicDataQuality.TPCFlag
        self._ECALFlag[0] = self.BasicDataQuality.ECALFlag
        self._DSECALFlag[0] = self.BasicDataQuality.DSECALFlag

        # Flag events from my bad_run_list
        if self._RunID[0] in self._bad_run_list:
            self.tianluFlag = True

        # Check to see whether file is run 1 or run 2 in order to use correct
        # bunch structure
        if not self.iscosmic:
            self.runNum = P0DCutUtils.getRunNum(self._IsMC[0], self._RunID[0])
        else:
            self.runNum = -1

        # First fill the truth info from the GRooTrackerVertices
        self.NPIDs = self.Global.NPIDs
        if self._IsMC[0]:
            self.NTraj = self.Trajectories.NTraj
            self.NVert = self.Vertices.NVtx

        # Sort pids in spill into bunches
        (self.maxMomIndices,
         self.maxMomNegIndices,
         self.maxMomPosIndices,
         self.pidsPerBunch,
         self.pidIdxsPerBunch,
         self.p0dPidsPerBunch,
         self.p0dPidIdxsPerBunch) = P0DCutUtils.sortToBunches(
             self.PIDs,
             self.NPIDs,
             self._IsMC[0],
             self._EventTime[0],
             self.runNum)

        # PoT tree variables fill
        self._SpillTime[0] = self.BasicHeader.EventTime
        # keep track of the number of global PIDs per bunch
        for i, pidsInBunchI in enumerate(self.pidsPerBunch):
            self._PIDsPerBunchArray[i] = pidsInBunchI
        # This keeps track of the number of bunches that has at least 1 PID for
        # each spill
        self._nBunchesWPIDPerSpill[0] = len(
            filter(lambda n: n > 0, self.pidsPerBunch))

        self._ResetFlags()


    def CheckBeamSummaryAddPOT(self):
        """ This function is used to check if we have good beam summary data and
        adds POT as necesary
        """
        # If we are looking at data files, need to check BeamSummaryData
        if not self._IsMC[0]:
            if not self.BeamSummaryData.BeamSummaryDataStatus:
                return
            if self.BasicDataQuality.ND280OffFlag != 0 or self.BasicDataQuality.P0DFlag != 0:
                return
            if not self.BeamData:
                return
            if self.BeamData.GoodSpillFlag == 0:
                return
            # If Data spill passes these cuts, we keep it
            self._POT[0] += self.BeamData.CT5ProtonsPerSpill
            self._PCT5perSpill[0] = self.BeamData.CT5ProtonsPerSpill
            return True
        else:
            # MC keeps track of POT file by file.  For GENIE MC, POT per file stored
            # in GRooTrackerVtx.Vtx.OrigTreePOT, we multiply by number of trees in
            # our TChain because Ntrees = NFiles.  We only need to check it once because
            # we're assuming the input files are all of same type (i.e. no mixing of
            # Data/MC or GENIE/NEUT files)
            if P0DCutUtils.isGENIE(self._RunID[0]):
                # Since the trees are all TChained together, and the POT corresponds
                # to PoT per tree/file, we need to ensure we're only adding to POT
                # count (self._POT[0])  when we encounter the next tree in the
                # TChain
                if self.NextPOTTreeNumber == self.GRooTrackerVtx.GetTreeNumber():
                    self._POT[0] += self.GRooTrackerVtx.Vtx.At(0).OrigTreePOT
                    # print 'Tree Number: ' + str(GRooTrackerVtx.GetTreeNumber())
                    # print 'Total POT Count: ' + str(self._POT[0])
                    self.NextPOTTreeNumber += 1
            else:
                if self.NextPOTTreeNumber == self.NRooTrackerVtx.GetTreeNumber():
                    self._POT[0] += self.NRooTrackerVtx.Vtx.At(0).OrigTreePOT
#                     print 'Tree Number: ' + str(GRooTrackerVtx.GetTreeNumber())
#                     print 'Total POT Count: ' + str(self._POT[0])
                    self.NextPOTTreeNumber += 1

            # No CT5ProtonsPerSpill for MC files
            self._PCT5perSpill[0] = 0
            return True


    def FillVarsFromTruth(self, vtx_selection, traj_selection):
        """ This function will loop over the true vertices in the oaAnal
        TTruthVerticeModule and select all the true events that pass the
        vtx_selection parameter.  It will fill the truth variables to the 
        truthTree.
        """
        for vtxN in range(self.NVert):
            trueVtx = self.Vertices.Vertices[vtxN]
            if not vtx_selection(trueVtx):
                # If this trueVtx fails the selection function
                continue

            # Get the trajectory for outgoing particle info
            trueTraj = self._SelectSingleTraj(trueVtx, traj_selection)
            if not trueTraj:
                # Even if the vertex passed, if no trajectory was matched
                # then we cannot detect it in the detector, which means
                # it really shouldn't count against our efficiency.
                continue

            self._NeutrinoPDG[0] = trueVtx.NeutrinoPDG
            self._TargetPDG[0] = trueVtx.TargetPDG

            self._Generator.Clear()
            self._Generator.Append(trueVtx.Generator)

            reaction = trueVtx.ReactionCode
            self._ReactionCode.Clear()
            self._ReactionCode.Append(reaction)
            self._ReactionConverted[0] = self.UnifiedReactionCode(reaction)

            pos = trueVtx.Position
            self._VertPos.SetPxPyPzE(pos.Px(),
                                       pos.Py(),
                                       pos.Pz(),
                                       pos.E())

            neutrinoMom = trueVtx.NeutrinoMomentum
            self._NeutrinoMomentum.SetPxPyPzE(neutrinoMom.Px(),
                                              neutrinoMom.Py(),
                                              neutrinoMom.Pz(),
                                              neutrinoMom.E())

            targetMom = trueVtx.TargetMomentum
            self._TargetMomentum.SetPxPyPzE(targetMom.Px(),
                                            targetMom.Py(),
                                            targetMom.Pz(),
                                            targetMom.E())

            trajPos = trueTraj.InitPosition
            trajMom = trueTraj.InitMomentum

            self._InitPos.SetXYZT(trajPos.X(), trajPos.Y(),
                               trajPos.Z(), trajPos.T())

            self._InitMom.SetPxPyPzE(trajMom.Px(), trajMom.Py(),
                                  trajMom.Pz(), trajMom.E())
            self._TrajPDG[0] = trueTraj.PDG

            self.truthTree.Fill()


    def FillVarsFromSelectedPID(self, selected):
        """ We pass this a selected GlobalPID and this will fill all the relevant info
        from that PID
        """
        # if there's a selected PID then set goodGlobalTrack to True
        self._goodGlobalTrack[0] = True

        # Unique trajectory ID associated with selected;
        # for matching of tracks between P0D, Global,
        # and Truth.
        trajID = selected.TrueParticle.ID
        vtxID = selected.TrueParticle.Vertex.ID

        self._Charge[0] = selected.Charge
        # Get Global info
        particleIDs = selected.ParticleIds
        pidWeights = selected.PIDWeights
        globalPos = selected.FrontPosition
        globalDir = selected.FrontDirection
        globalCone = selected.Cone
        self._nP0Ds[0] = selected.NP0Ds
        self._nTPCs[0] = selected.NTPCs
        self._nPDGs[0] = particleIDs.size()

        # keep track of which bunch the event fell
        self._BunchNum[0] = P0DCutUtils.findBunch(
            self._IsMC[0],
            self._EventTime[0],
            self.runNum,
            globalPos.T())
        # If a bunch was not found for this pid, the bunch num is set to -1.
        if self._BunchNum[0] != -1:
            # Get the number of GlobalPIDs in this event's bunch
            self._PIDsInBunch[0] = self.pidsPerBunch[self._BunchNum[0]]
            # Get the number of GlobalPIDs that have a p0d component
            self._P0DPIDsInBunch[0] = self.p0dPidsPerBunch[self._BunchNum[0]]
            # increment pot_tree counter that keeps track of n-selected PIDs
            # per bunch
            self._SelectedPIDsPerBunchArray[self._BunchNum[0]] += 1

        if self._nP0Ds[0] > 0:
            # Get P0DObject associated with globalPID
            # there should only be 1!
            p0dObj = selected.P0D.At(0)
            p0dPos = p0dObj.FrontPosition
            p0dBPos = p0dObj.BackPosition
            p0dFMom = p0dObj.FrontMomentum
            p0dBMom = p0dObj.BackMomentum
            p0dCone = p0dObj.Cone
            p0dFDir = p0dObj.FrontDirection

            # Fill Vars
            self._P0DLen[0] = p0dObj.Length
            self._P0DPosDiff[0] = (p0dPos - p0dBPos).Vect().Mag()
            self._P0DFPos.SetXYZT(p0dPos.X(), p0dPos.Y(),
                                  p0dPos.Z(), p0dPos.T())
            self._P0DBPos.SetXYZT(p0dBPos.X(), p0dBPos.Y(),
                                  p0dBPos.Z(), p0dBPos.T())
            self._P0DFMom[0] = p0dFMom
            self._P0DBMom[0] = p0dBMom
            self._P0DCone.SetXYZ(p0dCone.X(), p0dCone.Y(), p0dCone.Z())
            self._P0DFDir.SetXYZ(p0dFDir.X(), p0dFDir.Y(), p0dFDir.Z())

        # Now we check whether globalPID uses TPC1
        self._TPC1Used[0] = False
        for idx in range(self._nTPCs[0]):
            aUsedTPC = selected.TPC.At(idx)
            # error checking (weird shit be breaking randomly)
            try:
                aUsedTPC.Detector
            except AttributeError:
                print 'Base.py: WTFFF aUsedTPC is not actually TPC object'
                print aUsedTPC, self._RunID[0], self._EventID[0]
                continue

            if aUsedTPC.Detector == 1:
                self._TPC1Used[0] = True
                self._TPC1NNodes[0] = aUsedTPC.NNodes
                self._TPC1Len[0] = aUsedTPC.Length
                tpc1Pos = aUsedTPC.FrontPosition
                tpc1Dir = aUsedTPC.FrontDirection
                self._TPC1FPos.SetXYZT(tpc1Pos.X(), tpc1Pos.Y(),
                                       tpc1Pos.Z(), tpc1Pos.T())
                self._TPC1FMom[0] = aUsedTPC.FrontMomentum
                self._TPC1FDir.SetXYZ(tpc1Dir.X(), tpc1Dir.Y(),
                                      tpc1Dir.Z())
                self._TPC1Charge[0] = aUsedTPC.Charge
                self._TPC1Chi2[0] = aUsedTPC.Chi2
                break

        for j in range(particleIDs.size()):
            # pidPDGs.append((particleIDs[j], pidWeights[j]))
            self._PDGs[j] = int(particleIDs[j])
            self._PDGWeights[j] = float(pidWeights[j])
            # pidPDGs.append((particleIDs[j], pidWeights[j]))
            # print truePDG, pidPDGs
        self._GlobalFPos.SetXYZT(globalPos.X(), globalPos.Y(),
                                 globalPos.Z(), globalPos.T())
        self._GlobalFDir.SetXYZ(globalDir.X(), globalDir.Y(),
                                globalDir.Z())
        self._GlobalFMom[0] = selected.FrontMomentum
        self._GlobalCone.SetXYZ(globalCone.X(), globalCone.Y(),
                                globalCone.Z())
        self._isForward[0] = selected.isForward

        if self._IsMC[0]:
            # Get Truth info
            self._Eff[0] = selected.TrueParticle.Eff
            self._Pur[0] = selected.TrueParticle.Pur
            # Find matching truth track
            truthMatchBool = False
            for trajN in range(self.NTraj):
                trueTraj = self.Trajectories.Trajectories[trajN]
                if trueTraj.ID == trajID:
                    truthMatchBool = True
                    trajPos = trueTraj.InitPosition
                    trajMom = trueTraj.InitMomentum

                    self._FPos.SetXYZT(trajPos.X(), trajPos.Y(),
                                       trajPos.Z(), trajPos.T())

                    self._FMom.SetPxPyPzE(trajMom.Px(), trajMom.Py(),
                                          trajMom.Pz(), trajMom.E())

                    truePDG = trueTraj.PDG
                    self._TruePDG[0] = int(truePDG)
            if not truthMatchBool:
                print 'No true traj match!'

            # Match True vertex
            truthVtxBool = False
            for vtxN in range(self.NVert):
                trueVtx = self.Vertices.Vertices[vtxN]
                if trueVtx.ID == vtxID:
                    truthVtxBool = True

                    vtxPos = trueVtx.Position
                    # incoming neutrino info
                    neutrinoMom = trueVtx.NeutrinoMomentum
                    nuPDG = trueVtx.NeutrinoPDG
                    reaction = trueVtx.ReactionCode

                    self._VtxPos.SetXYZT(vtxPos.X(), vtxPos.Y(),
                                         vtxPos.Z(), vtxPos.T())
                    self._TrueVtxInP0DFiducial[0] = bool(
                        P0DCutUtils.inP0DFiducial(vtxPos))
                    self._NuMom.SetPxPyPzE(neutrinoMom.Px(), neutrinoMom.Py(),
                                           neutrinoMom.Pz(), neutrinoMom.E())
                    self._NuPDG[0] = nuPDG

                    # Use true nuPDG to determine FluxRatio for reweighting
                    # neutrino energy in GeV
                    nuEnergy = neutrinoMom.E() / 1000.
                    if nuPDG == 14:
                        self._FluxRatio[0] = (
                            self.fluxRatio_numu).Interpolate(nuEnergy)
                    elif nuPDG == -14:
                        self._FluxRatio[0] = (
                            self.fluxRatio_numub).Interpolate(nuEnergy)
                    elif nuPDG == 12:
                        self._FluxRatio[0] = (
                            self.fluxRatio_nue).Interpolate(nuEnergy)
                    elif nuPDG == -12:
                        self._FluxRatio[0] = (
                            self.fluxRatio_nueb).Interpolate(nuEnergy)
                    else:
                        self._FluxRatio[0] = 1

                    self._Reaction.Clear()
                    self._Reaction.Append(reaction)
                    self._ReactionInt[0] = self.UnifiedReactionCode(reaction)

            if not truthVtxBool:
                print 'No true vtx match!'
                self._FluxRatio[0] = 1

            self._trueGlobalMatched[0] = truthMatchBool
            self._vtxGlobalMatched[0] = truthVtxBool
        else:
            # No flux reweighting for Data
            self._FluxRatio[0] = 1


    def UnifiedReactionCode(self, reaction):
        """ cosmics, GENIE and NEUT use different reaction coding
        methods.  This function unifies them to a common integer form
        """
        if self.iscosmic:
            return int(reaction.split(' ')[1])
        elif P0DCutUtils.isGENIE(self._RunID[0]):
            # GENIE
            return P0DCutUtils.reactionStrToInt(reaction)
        else:
            # NEUT
            return P0DCutUtils.reactionCodeToInt(int(reaction))


    def MatchSelectedGlobal(self, selected_sub_obj, clonesarray):
        """ Pass me a selected globalPID subdetector object and I will match to clonesarray, if it exists, by matching
        unique IDs.  This is useful if we want to find the corresponding subdetector object that matches
        the selected globalpid.

        selected: TGlobalPID::TP0DObject (e.g. globalpid.P0D.At(0))
        clonesarray: TClonesArray of TP0DRecon:TP0DTrack/Shower objects (anything with a NHits specifier)
        """

        # first check that the objects in clonesarray have a UniqueID specifier
        try:
            global_unique = selected_sub_obj.NHits
            clonesarray.At(0).NHits
        except AttributeError:
            # print 'Base.py: Attribute error, cannot find NHits'
            return None
        except:
            print "Unexpected error:", sys.exc_info()[0]
            return None

        # try to match a p0dreco track to global_p0d_obj
        nEntries = clonesarray.GetEntriesFast()
        for i in range(nEntries):
            a_clones_obj = clonesarray.At(i)
            a_unique = a_clones_obj.NHits
            if a_unique == global_unique:
                return a_clones_obj

        return None


    def ProcessP0DRecoNodes(self, nodevec, init_mom, backwardsgoing, p0d2tpc):
        """ nodevec is a vector of integers (shorts) containing the indices
        to the nodes
        backwardsgoing: set to true if we believe the sample is a backwards going
        muon
        p0d2tpc: set to true if we want to start with initmom=0 at the front pos
        of the track
        """
        listofnodes = list(nodevec)
        curr_mom = init_mom
        if not p0d2tpc:
            # nodes are ordered from upstream->downstream.  This is good
            # for calculating the energy loss starting at 0, but we need a backwards
            # list of nodes if we want to start at TPC1
            listofnodes.reverse()

        # print 'listofnodes -------'
        # print listofnodes
        for i in range(len(listofnodes)):
            curr_node = self.P0D.Nodes.At(listofnodes[i])
            if i == len(listofnodes) - 1:
                # We are on the last node so break
                break

            next_node = self.P0D.Nodes.At(listofnodes[i + 1])

            # length between current and next node
            # strictly pos
            dl = (next_node.Position - curr_node.Position).Vect().Mag()

            # if we are looking at backwards going tracks and NOT using the
            # start with initmom=0 method, then path-length must be passed to
            # the ecorrector as negative
            if backwardsgoing and not p0d2tpc:
                dl = -dl

            if p0d2tpc:
                endz = next_node.Position.Z()
            else:
                endz = curr_node.Position.Z()

            curr_mom = self.eCorr.Calc(
                dl,
                curr_mom,
                curr_node.Direction,
                p0d2tpc,
                endz)

            # #### TESTING ####
            # print 'curr_node--next_node', curr_node.Position.Z(), next_node.Position.Z()
            # print 'initmom... currmom', init_mom, '...', curr_mom

        return curr_mom


    def InitHasPassedCutDict(self):
        """ initialize dictionary to keep track of whether an event has passed
        a registered cut
        """
        hasPassedCut = {}
        for step in self.listOfAnalysisSteps:
            hasPassedCut[step] = False

        hasPassedCut['No Cuts'] = True
        return hasPassedCut


    def MatchTPC2P0D(self):
        """ The purpose of this method is the match the TPC1 Object as given by
        the GlobalPID to a P0DRecon object.
        """
        self._P0DGoing[0] = self._CheckP0DGoing()

        # Only try to match TPC1 and P0DRecon if the track is seen to be
        # P0D-going
        if self._P0DGoing[0]:
            # Get the mcmtimesincebusy units of [10ns]
            self._MCMTimeSinceBusy[0] = self._GetMCMTimeSinceBusy()
            # return if the MCMTSB flag is -1 as that indicates the MCMTSB is
            # not in the dict (probably due to missing reco file)
            if self._MCMTimeSinceBusy[0] < 0:
                # reset flags as we basically want to ignore this event so it
                # shouldn't count as P0DGoing
                self._ResetFlags()
                return

            # This formula is based on the assumption that the MCMTSB has a zero
            # corresponding to the beginning of a series of P0D-cycles each
            # lasting 580 ns.  Since the cosmics trigger occurs randomly, we
            # need to adjust the timing relative to that to be relative to the
            # P0D cycles.
            # Two methods to calculate adjusted hit-times.
            f = lambda hit_time, mcmtsb: (hit_time + mcmtsb * 10) % 580
            f2 = lambda hit_time, mcmtsb, cap: (
                hit_time + mcmtsb * 10) % 580 + cap * 580

            # Get nP0DHits
            self._nP0DHits[0] = self.P0D.NHits
            # PyROOT doesn't seem to play nice with arrays that change size.
            # self._ResizeP0DHitArrsIfNecessary()

            n_good_hits = 0
            for i in range(self._nP0DHits[0]):
                p0dhit = self.P0D.Hits.At(i)
                self._P0DHitCharges[i] = p0dhit.Charge
                self._P0DHitTimes[i] = p0dhit.Time

                # Bitwise shift to get the capacitor (cycle) number.
                # Clark's suggestion
                self._P0DHitChanID[i] = p0dhit.ChanID
                cap = (self._P0DHitChanID[i] >> 7) & 0x1F

                # adjust hit time by offset
                # Some hit-times are set at ~4E7
                # Need to ensure that p0dhit.time is valid (i.e. has Good TDC)
                if p0dhit.Time < 10000:
                    adjusted_hit_time = f(p0dhit.Time,
                                          self._MCMTimeSinceBusy[0])
                    adjusted_hit_time_2 = f2(p0dhit.Time,
                                             self._MCMTimeSinceBusy[0],
                                             cap)
                else:
                    adjusted_hit_time = -p0dhit.Time
                    adjusted_hit_time_2 = -p0dhit.Time

                self._P0DAdjHitTimes[i] = adjusted_hit_time
                self._P0DAdjHitTimes2[i] = adjusted_hit_time_2

                n_good_hits += self._GoodHit(p0dhit, adjusted_hit_time)

            self._nGoodHits[0] = n_good_hits
            self._P0DActive[0] = n_good_hits >= 3

            # Get offset-adjusted TPC1FPos.T()
            self._TPC1AdjTime[0] = f(
                self._TPC1FPos.T(),
                self._MCMTimeSinceBusy[0])

            self._P0DReconTrackMatch[0] = self._P0DRecoTrackMatches()

            self._P0DReconShowerMatch[0] = self._P0DRecoShowerMatches()

            ## DEBUGGING ##
            ## missing p0d reco tracks/showers ##
            # if self._P0DActive[0] and not (self._P0DReconTrackMatch[0] or self._P0DReconShowerMatch[0]):
            #     for i in range(self.P0D.NHits):
            #         p0d_hit = self.P0D.Hits.At(i)
            # print p0d_hit.Time, p0d_hit.Charge, p0d_hit.ChanID,
            # p0d_hit.GeomID

            #     for i in range(self.P0D.NAlgoResults):
            #         algo_res = self.P0D.AlgoResults.At(i)
            #         '----'
            #         print algo_res.FullName
            #         print 'Cycle:', algo_res.Cycle
            ## END ##

            # Check if global matched a p0d track or shower
            if self._nP0Ds[0] > 0:
                self._P0DGlobalTrackMatch[0] = self._P0DCone.Mag() == 0
                self._P0DGlobalShowerMatch[0] = self._P0DCone.Mag() > 0


    ### Private methods below
    def _ResetFlags(self):
        """ Reset flags to default
        """
        self._trueGlobalMatched[0] = False
        self._vtxGlobalMatched[0] = False
        self._goodGlobalTrack[0] = False
        self._P0DGoing[0] = False
        self._P0DActive[0] = False
        self._nGoodHits[0] = 0
        self._nP0DHits[0] = 0
        self._P0DReconTrackMatch[0] = False
        self._P0DReconShowerMatch[0] = False
        self._P0DGlobalTrackMatch[0] = False
        self._P0DGlobalShowerMatch[0] = False
        self._MCMTimeSinceBusy[0] = 0


    def _SelectSingleTraj(self, trueVtx, selection):
        """ Returns the first true trajectory that is in the trueVtx's primary
        trajectories list AND passes the selection function.
        """

        n_vtx_traj = trueVtx.PrimaryTrajIDs.size()
        curr_vtxtraj_idx = 0
        # we can use the fact that the vtx primary IDs and the traj primary
        # IDs are sorted in ascending order to match trajectories using only
        # one for loop.  We need to keep track of when a match is found and
        # increment the curr_vtxtraj_idx accordingly.
        for atraj_idx in range(self.NTraj):
            atraj = self.Trajectories.Trajectories.At(atraj_idx)
            a_vtx_trajID = trueVtx.PrimaryTrajIDs.at(curr_vtxtraj_idx)

            ### DEBUGGING ###
            # print '--------------------------'
            # for i in range(trueVtx.PrimaryTrajIDs.size()):
            #     print ('traj id, vtx_traj_id', atraj.ID,
            #            atraj.PrimaryID, trueVtx.PrimaryTrajIDs.at(i),
            #            atraj.PDG)
            ###### END ######

            # While atraj.PrimaryID is larger than the current vtx traj ID, we
            # need to select the next trajectory associated with the vtx as it
            # will have a larger PrimaryID.  Only when the vtx_traj's Primary
            # ID is larger than or equal to atraj.PrimaryID can we start
            # checking whether the trajectory IDs match or continue to
            # increment atraj_idx
            while (atraj.PrimaryID > a_vtx_trajID):
                if curr_vtxtraj_idx+1 >= n_vtx_traj:
                    # this means we have run out of trajectories associated
                    # with the vtx to check.
                    print 'Base.py: No Match!', self._RunID[0], self._SubrunID[0], self._EventID[0]
                    return None

                curr_vtxtraj_idx += 1
                a_vtx_trajID = trueVtx.PrimaryTrajIDs.at(curr_vtxtraj_idx)

            if atraj.PrimaryID < a_vtx_trajID:
                # if the trajectory's primaryID is less than the current
                # vtx trajectory's primaryID, there is no match and we go to
                # the next trajectory
                continue
            elif (atraj.PrimaryID == a_vtx_trajID
                  and selection(atraj)):
                # if the PrimaryIDs match, two things can happen.  Either
                # the trajectory passes the selection, in which case it is
                # returned as the Selected Trajectory, or it does not pass
                # the selection and nothing is returned for now, but we must
                # increment the curr_vtxtraj_idx as the one associated with
                # atraj doesn't pass our check.
                return atraj
            elif curr_vtxtraj_idx+1 < n_vtx_traj:
                # at this point, we know that
                # atraj.PrimaryID == a_vtx_trajID but 
                # selection(atraj) returned False.  Need to increment if
                # we still have trajectories associated with the vertex.
                curr_vtxtraj_idx += 1
            else:
                # At this point there are no trajectories associated with
                # the vertex left to check
                print 'Base.py: No match 2!', self._RunID[0], self._SubrunID[0], self._EventID[0]
                return None

        # at this point there are no trajectories left
        print 'Base.py: No traj left', self._RunID[0], self._SubrunID[0], self._EventID[0]
        return None


    ### Functions for matching TPC objects to P0DRecon and generating mcmtsbdict
    def _CreateMarshaledFile(self, **kwargs):
        msh_file = kwargs['marshal_file']
        lock = FileLock(msh_file)
        lock.acquire()
        ch = ROOT.TChain(kwargs['tree_name'])
        for chain_file in kwargs['to_chain_files']:
            ch.AddFile(chain_file)

        vars_to_get = kwargs['var_to_get']

        self._mcm_time_dict = P0DCutUtils.convertToDict(ch, **vars_to_get)

        with open(msh_file, 'wb') as f:
            marshal.dump(self._mcm_time_dict, f)

        lock.release()


    def _SetupMCMTimeDictionary(self):
        to_chain_base = directories.TIANND280_STORE
        kwargs = {'var_to_get': {'MCMTimeSinceBusy': 'i'},
                  'tree_name': 'getmcmtimesTree',
                  'to_chain': os.path.join(to_chain_base,
                                           'getMCMTimes/*/*/*/*/*.root'),
                  'marshal_file': os.path.join(to_chain_base,
                                               'getMCMTimes/mcm_dict.out')}

        marshal_file = kwargs['marshal_file']
        to_chain_files = glob.glob(kwargs['to_chain'])
        kwargs['to_chain_files'] = to_chain_files

        if not to_chain_files:
            print 'Basy.py: No files to chain for MCMTime!'
            return

        chain_mtimes = [os.path.getmtime(f) for f in to_chain_files]

        # the marshal_file is locked only when CreateMarshaledFile is run.
        # If it is, then some process is runnign CreateMarshaledFile(), so
        # we need to wait until that process finishes before continuing.
        lock = FileLock(marshal_file)
        while lock.is_locked():
            time.sleep(10)

        # if the marshal_file (binary where dict is saved) doesn't exist,
        # create it
        if not os.path.exists(marshal_file):
            print 'Creating marshaled file...'
            self._CreateMarshaledFile(**kwargs)
        # Else if the marshal_file is older than one of the *.root files
        # create it
        elif os.path.getmtime(marshal_file) < max(chain_mtimes):
            print 'Updating marshaled file...'
            self._CreateMarshaledFile(**kwargs)
        # Otherwise, load up marshaled file and save to self._mcm_time_dict
        else:
            print 'Loading marshaled file...'
            with open(marshal_file, 'rb') as msh_file:
                self._mcm_time_dict = marshal.load(msh_file)


    def _GetMCMTimeSinceBusy(self):
        """ Get the MCMTimeSinceBusy, which is used to calculate the timing edges for
        cosmics and such.
        """
        # Only need the mcmtime for data cosmics
        # return 0 if the event is a MC event or is not a cosmic.
        # Will use MCMTSB==0 to indicate all the hits have valid times
        if self._IsMC[0] or not self.iscosmic:
            return 0

        if not self._mcm_time_dict:
            self._SetupMCMTimeDictionary()

        try:
            key_base = '{run}\{event}\{var}'
            key_name = key_base.format(run = self._RunID[0],
                                       event = self._EventID[0],
                                       var = 'MCMTimeSinceBusy')

            return self._mcm_time_dict[key_name]
        except KeyError:
            print 'Base.py:', self._RunID[0], self._EventID[0], 'has no corresponding mcmtime'
            return -1


    def _CheckP0DGoing(self):
        """ Check to see if the event is a P0DGoing candidate
        """

        # Require TPC1 be used
        if not self._TPC1Used[0]:
            return False
        # Require that TPC1FPos.Z() < -755
        if self._TPC1FPos.Z() >= -755:
            return False
        # Require > 18 nodes.  This is the TPC track quality cut, it rejects short tracks
        # for which the reconstruction is less reliable
        if self._TPC1NNodes[0] <= 18:
            return False

        # Project the track from the given TPC1FPosition to the P0D at a preset z-plane of -1100
        # Ensure that the projected position is within an XY fiducial volume
        proj_pos = P0DCutUtils.projectIntoP0D(self._TPC1FPos, self._TPC1FDir)
        in_xy_fid = (P0DCutUtils.XCUT_FIDUCIAL_MIN
                     < proj_pos.X() < P0DCutUtils.XCUT_FIDUCIAL_MAX and
                     P0DCutUtils.YCUT_FIDUCIAL_MIN
                     < proj_pos.Y() < P0DCutUtils.YCUT_FIDUCIAL_MAX)

        # if we're looking at cosmic triggers we require that they be
        # FGDCosmics
        if self.iscosmic:
            p0d_going = (self._FGDCosmic[0] and
                         self._TPC1FMom[0] > P0DCutUtils.TPC1FMOM_MIN and
                         in_xy_fid)
        # otherwise must be beam triggers, don't require cosmics trigger flag
        else:
            p0d_going = self._TPC1FMom[0] > P0DCutUtils.TPC1FMOM_MIN and in_xy_fid

        return p0d_going


    def _ValidHitTime(self, t):
        if self._MCMTimeSinceBusy[0] == 0:
            # if mc or not fgdcosmic, MCMTimeSinceBusy is manually set to 0
            # hit time is always valid
            return True
        elif self._MCMTimeSinceBusy[0] == -1:
            # if there was no corresponding mcmTimeSinceBusy, treat every
            # hit-time as invalid
            return False
        else:
            # Require hits be "edge" ns from the edges of a P0D
            # integration window
            # edge = 80
            # return -4880 + edge < t < -4400 - edge
            edge = 0
            return 0 + edge < t < 200 - edge


    def _GoodHit(self, p0dhit, adjusted_hit_time):
        # TODO: Probably should implement some sort of hit position check
        return self._ValidHitTime(adjusted_hit_time)


    def _MatchesTPC1FDir(self, p0d_dir, dir_stddev):
        upper = p0d_dir + dir_stddev + dir_stddev
        lower = p0d_dir - dir_stddev - dir_stddev

        return (lower.X() < self._TPC1FDir.X() < upper.X() and
                lower.Y() < self._TPC1FDir.Y() < upper.Y() and
                lower.Z() < self._TPC1FDir.Z() < upper.Z())


    def _ValidP0DRecoDir(self, a_reco_obj):
        """ Check that the direction associated with a P0DRecon Shower
        matches the TPC1FDir within 2-stddev.

        The direction can be either as reconstructed, or reflected about
        its origin.  I.e. negative.
        """
        dir_stddev = ROOT.TVector3(math.sqrt(a_reco_obj.DirVariance.X()),
                                   math.sqrt(a_reco_obj.DirVariance.Y()),
                                   math.sqrt(a_reco_obj.DirVariance.Z()))

        return (self._MatchesTPC1FDir(a_reco_obj.Direction, dir_stddev) or
                self._MatchesTPC1FDir(-a_reco_obj.Direction, dir_stddev))


    def _TPC1FDirInCone(self, a_reco_obj):
        """ Check that the angle between TPC1FDir and P0DRecon shower direction is less than the
        opening angle of the P0DShower::Cone.  Assume that shower direction can be mis-recon to
        be negative so we check both the inner-angle, and pi-inner-angle
        """

        # print '---'
        # print a_reco_obj
        # print a_reco_obj.Cone,
        # rootutils.inner_angle(a_reco_obj.Direction, self._TPC1FDir)
        cone_angle = a_reco_obj.Cone
        inner_angle = rootutils.inner_angle(a_reco_obj.Direction,
                                              self._TPC1FDir)
        return inner_angle < cone_angle or (math.pi - inner_angle) < cone_angle


    def _P0DRecoTrackMatches(self):
        # Check if P0DRecon reconstructed a track
        # Loop through tracks, see if last node is exiting
        for i in range(self.P0D.NTracks):
            a_reco_track = self.P0D.Tracks.At(i)

            # most downstream node of track
            last_node_idx = list(a_reco_track.Nodes)[-1]
            last_node = self.P0D.Nodes.At(last_node_idx)

            # From TN80
            if last_node.Position.Z() > -1016:
                return True

        return False


    def _P0DRecoShowerMatches(self):
        for i in range(self.P0D.NShowers):
            a_reco_shower = self.P0D.Showers.At(i)

            if self._TPC1FDirInCone(a_reco_shower):
                return True

        return False


    def _ResizeP0DHitArrsIfNecessary(self):
        # Check to see if array size for storing hit-info is large enough, if
        # not, resize
        while self._nP0DHits[0] > len(self._P0DHitChanID):
            self._P0DHitChanID.extend([0] * len(self._P0DHitChanID))
            self._P0DHitTimes.extend([0.0] * len(self._P0DHitTimes))
            self._P0DAdjHitTimes.extend([0.0] * len(self._P0DAdjHitTimes))
            self._P0DAdjHitTimes2.extend([0.0] * len(self._P0DAdjHitTimes2))
            self._P0DHitCharges.extend([0.0] * len(self._P0DHitCharges))
