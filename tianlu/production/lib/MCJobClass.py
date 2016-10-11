"""
Class for storing all settings required to generate the config file for production.
Bases everything off the pwd and staging/runInfo.card.
"""
import os
import sys
import string
from glob import glob
from tianlu import utils


class MCRunInfoClass:
    def __init__(self):
        self.setDefaults()
        self.setDetails()
        self.MissingConfigs = []


    def fillFromCard(self,cardname):
        input = open(cardname, 'r')
        lines = []
        dirDict = {}
        for line in input:
            dirDict[line.split(':')[0]] = line.split(':')[1].strip()

        keyCheck = ['store', 'input', 'output']
        dirKeys = dirDict.keys()
        for key in keyCheck:
            if key not in dirKeys:
                print 'Check', cardname, 'for proper format, could not find \'%s\' key'%key
                sys.exit()
        self.store = dirDict['store']
        self.input = dirDict['input']
        self.output = dirDict['output']
        self.startModule = 'nd280MC'
        self.passThruDirFromCard = ''
        self.numcDirFromCard = ''
        self.fluxDirFromCard = ''
        if 'startmodule' in dirKeys:
            self.startModule = dirDict['startmodule']
        if 'passthrudir' in dirKeys:
            self.passThruDirFromCard = os.path.expandvars(
                os.path.expanduser(dirDict['passthrudir']))
        if 'numcthrudir' in dirKeys:
            self.numcDirFromCard = os.path.expandvars(
                os.path.expanduser(dirDict['numcthrudir']))
        if 'fluxdir' in dirKeys:
            self.passThruDirFromCard = os.path.expandvars(
                os.path.expanduser(dirDict['fluxdir']))


    def setDetails(self):
        """ Order of setting these parameters matters
        """
        pwd = utils.get_cwd()
        production_idx = pwd.find('production')+13
        self.prodDir = pwd[0:production_idx]
        self.nd280Version = os.path.basename(os.getenv('ND280ROOT'))
        self.production = self.prodDir.split('production')[-1].strip('0')
        self.respin = utils.split_path(pwd[production_idx+1:])[0]
        self.nuMCType = utils.split_path(pwd[production_idx+1:])[2]
        self.fillFromCard('runInfo.card')

        self.usingNUCP = False
        if 'beam/' in pwd:
            self.beam = 'beam'
        if 'nue/' in pwd:
            self.beam = 'nue'
            self.nuType = 'nue'
        if 'run1/' in pwd:
            self.beam = 'run1'
            self.ecalPeriods = '1-2'
        if 'run2/' in pwd:
            self.beam = 'run2'
            self.tpcPeriods = 'runs2-3'
            self.ecalPeriods = '1-2'
        if 'run3/' in pwd:
            self.beam = 'run3'
            self.tpcPeriods = 'runs2-3'
            self.ecalPeriods = '3-4'
        if 'run4/' in pwd:
            self.beam = 'run4'
            self.tpcPeriods = 'runs2-3-4'
            self.ecalPeriods = '3-4'
        if 'run5/' in pwd:
            self.beam = 'run5'
            self.tpcPeriods = 'runs2-3-4'
            self.ecalPeriods = '3-4'
        if 'ccpiplus/' in pwd:
            self.beam = 'ccpiplus'
            self.nMesons = 0
            self.nLeptons = 1
            self.nMuMinus = 1
            self.nPiZero = 0
            self.nPiPlus = 1
            self.usingNUCP = True
        if 'ccpizero/' in pwd:
            self.beam = 'ccpizero'
            self.nMesons = 0
            self.nLeptons = 1
            self.nMuMinus = 1
            self.nPiZero = 1
            self.nPiPlus = 0
            self.usingNUCP = True
        if 'ncpiplus/' in pwd:
            self.beam = 'ncpiplus'
            self.nMesons = 0
            self.nLeptons = 0
            self.nMuMinus = 0
            self.nPiZero = 0
            self.nPiPlus = 1
            self.usingNUCP = True
        if 'ncpizero/' in pwd:
            self.beam = 'ncpizero'
            self.nMesons = 0
            self.nLeptons = 0
            self.nMuMinus = 0
            self.nPiZero = 1
            self.nPiPlus = 0
            self.usingNUCP = True
        if 'tpcgas/' in pwd:
            self.beam = 'tpcgas'
        if 'verify/' in pwd:
            self.verify = True
        if self.nuMCType == 'anti-genie':
            self.runprefix -= 10000000
        if 'genie' in pwd:
            self.mc = 'Genie'
            self.runprefix += 1000000
        self.respin = pwd[pwd.find(self.prodDir)+len(self.prodDir)+1:][0]
        if self.respin not in string.uppercase:
            print 'Respin', self.respin, 'doesn\'t appear to be an UPPER CASE LETTER'
        if '2010-11' in pwd:
            self.baseline = '2010-11'

        if 'magnet/' in pwd:
            self.runN = int(pwd[pwd.find('/run')+4])
            self.runprefix += (self.runN-1)*100000

        if 'water' in pwd:
            self.fill = 'water'
            self.p0dwater = 1
            self.runprefix += 10000
        if 'basket/' in pwd:
            self.fluxVolume = 'basket'
            self.fluxMasterVolume = 'Basket'
            self.fluxName = 'basket'
            self.runN = 2
            self.runprefix += 101000
            if 'nue/' in pwd:
                self.fluxName = 'Nue'
                self.runprefix += 1000
            elif 'ncpizero/' in pwd:
                self.fluxName = 'NC1pi0'
                self.runprefix += 2000
            elif 'ccpizero/' in pwd:
                self.fluxName = 'CC1pi0'
                self.runprefix += 3000
            elif 'ncpiplus/' in pwd:
                self.fluxName = 'NC1pi+'
                self.runprefix += 4000
            elif 'ccpiplus/' in pwd:
                self.fluxName = 'CC1pi+'
                self.runprefix += 5000
            elif 'ncpizerofgd/' in pwd:
                self.fluxName = 'NCpi0FGD'
                self.fluxMasterVolume = 'FGD1'
                self.runprefix += 6000
            elif 'ccpicoh/' in pwd:
                self.fluxName = 'CCpicoh'
                self.fluxMasterVolume = 'FGD1'
                self.runprefix += 7000
            elif 'tpcgas/' in pwd:
                self.fluxName = 'TPCGas'
                # set this to mask ND280 geometry
                # the self.standalone option can be set to a single ND280 detector
                # and overrides the baseline setting.  However, turns out that
                # setting master_volume to Half produces events only on argon so
                # we are using that instead.
                # self.standalone = 'TPC'
                self.fluxMasterVolume = 'Half'
                self.forceVolume = 'true'
                self.runprefix += 6000

        self.setBasePath()
        self.setNumcDir()
        self.setPassThruDir()
        self.setFluxDir()
        self.setFluxInfo()


    def setDefaults(self):
        self.beam = 'beam'
        self.nuType = 'beam'
        self.fill = 'air'
        self.p0dwater = 0
        self.production = 5
        self.startModule = 'nd280MC'
        self.respin = 'A'
        self.baseline = '2010-02'
        self.runN = 1
        self.fluxVolume = 'magnet'
        self.fluxMasterVolume = 'Magnet'
        self.fluxName = 'magnet'
        self.mc = 'Neut'
        self.runprefix = 90100000
        self.verify = False
        self.setCherryPickingDefaults()
        self.tpcPeriods = None
        self.ecalPeriods = None
        self.standalone = None
        self.forceVolume = 'false'


    def setCherryPickingDefaults(self):
        self.inputFileList = []
        self.nMesons = 0
        self.nLeptons = 0
        self.nMuMinus = 0
        self.nPiZero = 0
        self.nPiPlus = 0


    def setBasePath(self):
        basePath = utils.get_cwd()
        basePath = basePath[0:basePath.rfind('staging')-1]
        self.basePath = basePath


    def setNumcDir(self):
        if self.numcDirFromCard:
            # This should be set if we want to use a reco/cali/etc file
            # as input but still need the numc file for the oaAnal passthrough
            self.NumcDir = self.numcDirFromCard
        elif self.passThruDirFromCard:
            # If the numcDirFromCard variable is not set, assume that
            # passthrudirfromcard contains the numc dir
            self.NumcDir = self.passThruDirFromCard
        else:
            ptdir = self.basePath+'/numc'
            # just swap the respin in the path
            ndex = ptdir.find('/mcp')-1
            dir = ptdir[0:ndex]+'A'+ptdir[ndex+1:]
            self.NumcDir = dir


    def setPassThruDir(self):
        # first if it's set in the card, then use that
        if self.passThruDirFromCard:
            self.PassThruDir = self.passThruDirFromCard
        else:
            if self.respin == 'A':
                dir = self.basePath+'/numc'
            else:
                dir = self.NumcDir
            if self.usingNUCP and self.output == 'nd280':
                dir = self.basePath+'/nucp'

            # T.Y. 11/14/13.    Commented out to try to keep verify sample
            # completely self contained through all MC modules.  The
            # following code, if left uncommented, will search for GENIE
            # vector files in the production dir (non-verify sample).
            #
            # if 'verify' in dir:
            #       version = self.nd280Version
            #       start = dir[0:dir.find('verify')]
            #       end = dir[dir.find(version)+len(version)+1:]
            #       dir = start+end
            self.PassThruDir = dir


    def setFluxDir(self):
        if self.fluxDirFromCard:
            self.fluxDir = self.fluxDirFromCard
        elif os.path.isdir(self.prodDir+'/flux'):
            self.fluxDir = self.prodDir+'/flux'
        elif self.output == 'numc':
            print 'Flux directory not set properly for numc!'
            sys.exit(1)
        else:
            self.fluxDir = None


    def setFluxInfo(self):
        self.fluxTree = 'h3002'
        self.XSFile = 'gxspl-t2k-v2.8.0.xml'
        self.randomStart = 1
        if self.fluxVolume == 'basket':
            if self.nuMCType == 'anti-genie':
                self.fluxFile = os.path.basename(
                    utils.glob_newest(os.path.join(self.fluxDir,
                                                   'nu.nd5m_flukain.0-[0-9]*.root')))
            else:
                self.fluxFile = os.path.basename(
                    utils.glob_newest(os.path.join(self.fluxDir,
                                                   'nu.nd5_flukain.0-[0-9]*.root')))
            # Note different PoT/file for different basket samples
            # (T.Y. 6/16/2014)
            if self.fluxName == 'NC1pi+':
                self.fluxPOT = '5e19'
            elif self.fluxName == 'CC1pi+':
                self.fluxPOT = '1e19'
            elif self.fluxName == 'basket':
                self.fluxPOT = '1e18'
            elif self.fluxName == 'TPCGas':
                self.fluxPOT = '5e20'
            else:
                self.fluxPOT = '2e19'
        else:
            # magnet
            self.fluxPOT = '5e17'
            if self.beam == 'run5':
                # anti-nu flux
                self.fluxFile = os.path.basename(
                    utils.glob_newest(os.path.join(self.fluxDir,
                                                   'nu.nd6m_flukain.0-[0-9]*.root')))
            else:
                # nu
                self.fluxFile = os.path.basename(
                    utils.glob_newest(os.path.join(self.fluxDir,
                                                   'nu.nd6_flukain.0-[0-9]*.root')))

        geometry = (self.standalone.lower() if self.standalone
                    else self.baseline+'-'+self.fill)
        flux_type = os.path.splitext(self.fluxFile)[0]
        self.fluxProbs = flux_type+'_genie_evtrate_'
        self.fluxProbs += geometry+'_'+self.fluxMasterVolume.lower()+'.root'
        self.fluxFilePOT = str(int(flux_type.split('-')[-1])+1)+'e21'


    def getNeededFiles(self):
        if not self.XSFile:
            self.setFluxInfo()

        flux_probs_dir = 'genie_evtrate/'
        files = [self.XSFile, self.fluxFile, flux_probs_dir+self.fluxProbs]
        return files


    def getJobName(self, run, subrun):
        handle = 'R'+str(self.runN)+self.fluxMasterVolume[0].lower()+self.fill[0]
        handle += '_'+str(run)+'_'+str(subrun)
        return handle

