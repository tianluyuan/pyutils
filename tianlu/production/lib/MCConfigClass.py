#!/usr/bin/env python

import os
import sys
import getpass
import random
import glob
from tianlu.utils import queue_check, make_dirs_if_needed, get_cwd
import tianlu.directories as tdir

# system-dependent (filepath) stuff:
genieSetupPath = '/nfs/hepcode/t2k/GENIE/setup.sh'

#------------------------------------------------------------------------------
# MCConfigClass
#------------------------------------------------------------------------------

class MCConfigClass:
    def __init__(self,runInfo = None):
        if runInfo is None:
            print 'Please provide a valid RunInfo class'
        self.runInfo = runInfo
        self.cfgDict = {}
        self.setRunInfo()
        self.makeDicts()
        self.setSoftwareFields()
        self.setConfigFields()
        self.setFileNamingFields()
        self.setGeometryFields()
        self.setND280Fields()
        self.setCherryPickingFields()
        if runInfo.output == 'numc':
            #FIXME
            self.setNeutrinoFields()
        if runInfo.output == 'nd280':
            self.setSoftwareFields()
            self.setElectronicsCfgFields()
            self.setAnalysisFields()
            self.setDeadChannelFields()
            self.throwAgain()
        if runInfo.output == 'nucp':
            self.cfgDict['module_list'] = 'oaCherryPicker'
            self.cfgDict['num_events'] =    100000000
            self.throwAgain()


    def throwAgain(self):
        self.cfgDict['mc_seed'] = random.randint(0,9e8)
        self.cfgDict['elec_seed'] = random.randint(0,9e8)
        self.mcSeed = random.randint(0,9e8)
        self.elecSeed = random.randint(0,9e8)



#------------------------------------------------------------------------------
# Initial Setters
#------------------------------------------------------------------------------

    def setRunInfo(self):
        self.cfgDict['p0d_water_fill'] = 1
        self.mcPosition = 'free'
        if self.runInfo.fill == 'air':
            self.cfgDict['p0d_water_fill'] = 0
        self.cfgDict['mc_type'] = self.runInfo.mc
        if self.runInfo.fluxVolume == 'basket':
            self.nBunches = 1
            self.bunchDuration = 0
            self.timeOffset = 0
            self.countType = 'FIXED'
            self.interactionsPerSpill = 1
            self.potPerSpill = 0
            self.mcFullSpill = 0

        if self.runInfo.fluxVolume == 'magnet':
            # NOTE: MUST CHANGE INTERACTIONSPERSPILL WHEN GENIE FLUXPROBS
            # OR POT/Sp VALUES CHANGE (T.Y)
            # Current values for prod6 using GENIE 2.8.0.  Should correspond
            # to values on mcProdSummary page
            # TODO: make this dynamically read from the interaction logs...smfh
            self.timeOffset = 50
            self.countType = 'MEAN'
            self.mcFullSpill = 1
            if self.runInfo.runN == 1:
                self.nBunches = 6
                self.potPerSpill = 3.66e13
                self.bunchDuration = 17
                self.interactionsPerSpill = 3.72428729
            if self.runInfo.runN == 2:
                self.nBunches = 8
                self.potPerSpill = 7.99e13
                self.bunchDuration = 19
                if self.runInfo.fill == 'water':
                    self.interactionsPerSpill = 8.51893711
                if self.runInfo.fill == 'air':
                    self.interactionsPerSpill = 8.49616818
            if self.runInfo.runN == 3:
                self.nBunches = 8
                self.potPerSpill = 9.462526e13
                self.interactionsPerSpill = 10.0636166
                self.bunchDuration = 19
            if self.runInfo.runN == 4:
                self.nBunches = 8
                self.potPerSpill = 9.462526e13
                self.bunchDuration = 19
                if self.runInfo.fill == 'water':
                    self.interactionsPerSpill = 10.0905862
                if self.runInfo.fill == 'air':
                    self.interactionsPerSpill = 10.0636166
            if self.runInfo.runN == 5:
                self.nBunches = 8
                self.potPerSpill = 11.3303e13
                self.bunchDuration = 19
                self.interactionsPerSpill = 4.26512148


    def createDirDict(self, dir_path):
        dir_dict = {}
        dir_files = glob.glob(dir_path+'/*.root')
        for aFile in dir_files:
            # key is the run-subrun identifier
            key = aFile.split('_')[3]
            dir_dict[key] = aFile

        return dir_dict


    def makeDicts(self):
        # create the dicts for numc and passthru dirs
        self.NumcDirFiles = self.createDirDict(self.runInfo.NumcDir)
        if self.runInfo.NumcDir != self.runInfo.PassThruDir:
            self.PassThruDirFiles = self.createDirDict(self.runInfo.PassThruDir)
        else:
            self.PassThruDirFiles = self.NumcDirFiles
#------------------------------------------------------------------------------
# Set the nd280Control fields
#------------------------------------------------------------------------------

    #-----------------
    # [software]
    #-----------------
    def setSoftwareFields(self):
        self.cfgDict['[software]'] = []
        self.cfgDict['[software]'].append(('cmtpath', 'environment'))
        self.cfgDict['[software]'].append(('cmtroot', 'environment'))
        self.cfgDict['[software]'].append(('nd280ver', self.runInfo.nd280Version))
        if self.runInfo.output == 'numc':
            self.cfgDict['[software]'].append(('genie_setup_script', genieSetupPath))

    #-----------------
    # [configuration]
    #-----------------
    def setConfigFields(self):
        self.inputFile = None
        self.cfgDict['[configuration]'] = []
        modList = ''
        if self.runInfo.output == 'numc':
            modList = 'genieMC genieConvert'
        if self.runInfo.output == 'nd280':
            modList = ''
            mods = ['nd280MC', 'elecSim', 'oaCalibMC', 'oaRecon', 'oaAnalysis']
            start = mods.index(self.runInfo.startModule)
            mods = mods[start:]
            for mod in mods:
                modList += mod+' '
        if self.runInfo.output == 'nucp':
            modList = 'oaCherryPicker'
        self.cfgDict['[configuration]'].append(('module_list', modList))

    #-----------------
    # [filenaming]
    #-----------------
    def setFileNamingFields(self):
        self.cfgDict['[filenaming]'] = []
        self.cfgDict['[filenaming]'].append(('run_number', None))
        self.cfgDict['[filenaming]'].append(('subrun', None))
        self.setComment()
        self.cfgDict['[filenaming]'].append(('comment', self.runInfo.Comment))

    #-----------------
    # [neutrino]
    #-----------------
    def setNeutrinoFields(self):
        self.NeutrinoCfg = True
        self.cfgDict['[neutrino]'] = []
        fluxVolume = self.runInfo.fluxVolume
        fluxMasterVolume = self.runInfo.fluxMasterVolume
        nuType = self.runInfo.nuType
        self.cfgDict['[neutrino]'].append(('neutrino_type',nuType))
        self.cfgDict['[neutrino]'].append(('master_volume', fluxMasterVolume))
        self.cfgDict['[neutrino]'].append(('force_volume_name', self.runInfo.forceVolume))
        self.cfgDict['[neutrino]'].append(('flux_region', fluxVolume))
        self.cfgDict['[neutrino]'].append(('flux_file', self.runInfo.fluxFile))
        self.cfgDict['[neutrino]'].append(('flux_tree', self.runInfo.fluxTree))
        self.cfgDict['[neutrino]'].append(('flux_file_pot', self.runInfo.fluxFilePOT))
        self.cfgDict['[neutrino]'].append(('pot', self.runInfo.fluxPOT))
        self.cfgDict['[neutrino]'].append(('genie_xs_table', self.runInfo.XSFile))
        self.cfgDict['[neutrino]'].append(('random_start', self.runInfo.randomStart))
        self.cfgDict['[neutrino]'].append(('genie_flux_probs_file_name',
                                           self.runInfo.fluxProbs))
    #-----------------
    # [analysis]
    #-----------------
    def setAnalysisFields(self):
        save_geom_opt = str(int('6' == self.runInfo.production))
        self.cfgDict['[analysis]'] = []
        self.cfgDict['[analysis]'].append(('save_geometry', save_geom_opt))
        #self.cfgDict['[analysis]'].append(('pass_through_dir', self.runInfo.PassThruDir))
        #self.cfgDict['NeutrinoCfg'].append(('random_seed',

    #-----------------
    # [geometry]
    #-----------------
    def setGeometryFields(self):
        self.cfgDict['[geometry]'] = []
        self.cfgDict['[geometry]'].append(('baseline', self.runInfo.baseline))
        self.cfgDict['[geometry]'].append(('p0d_water_fill', self.runInfo.p0dwater))
        if self.runInfo.standalone:
            self.cfgDict['[geometry]'].append(('standalone', self.runInfo.standalone))

    #-----------------
    # [electronics]
    #-----------------
    def setElectronicsCfgFields(self):
        self.cfgDict['[electronics]'] = []

    #-----------------
    # [nd280mc]
    #-----------------
    def setND280Fields(self):
        self.cfgDict['[nd280mc]'] = []
        self.cfgDict['[nd280mc]'].append(('mc_type', 'Genie'))
        if self.runInfo.output in ['nd280','nucp']:
            self.cfgDict['[nd280mc]'].append(('num_events', '100000000'))
            self.cfgDict['[nd280mc]'].append(('mc_full_spill', self.mcFullSpill))
            self.cfgDict['[nd280mc]'].append(('interactions_per_spill',
                                              self.interactionsPerSpill))
            self.cfgDict['[nd280mc]'].append(('time_offset', self.timeOffset))
            self.cfgDict['[nd280mc]'].append(('count_type', self.countType))
            self.cfgDict['[nd280mc]'].append(('mc_position', self.mcPosition))

        if self.runInfo.output == 'nd280':
            self.cfgDict['[nd280mc]'].append(('nbunches', self.nBunches))
            self.cfgDict['[nd280mc]'].append(('pot_per_spill', self.potPerSpill))
            self.cfgDict['[nd280mc]'].append(('bunch_duration', self.bunchDuration))


    #-----------------
    # [cherry_picker]
    #-----------------
    def setCherryPickingFields(self):
        self.cfgDict['[cherry_picker]'] = []
        if self.runInfo.output == 'nucp':
            self.cfgDict['[cherry_picker]'].append(('num_mesons', self.runInfo.nMesons))
            self.cfgDict['[cherry_picker]'].append(('num_leptons', self.runInfo.nLeptons))
            self.cfgDict['[cherry_picker]'].append(('num_mu_minus', self.runInfo.nMuMinus))
            self.cfgDict['[cherry_picker]'].append(('num_pizero', self.runInfo.nPiZero))
            self.cfgDict['[cherry_picker]'].append(('num_piplus', self.runInfo.nPiPlus))


    #-----------------
    # [dead_channels]
    #-----------------
    def setDeadChannelFields(self):
        self.cfgDict['[dead_channels]'] = []
        if self.runInfo.tpcPeriods:
            self.cfgDict['[dead_channels]'].append(('tpc_periods_to_activate', self.runInfo.tpcPeriods))
        if self.runInfo.ecalPeriods:
            self.cfgDict['[dead_channels]'].append(('ecal_periods_to_activate', self.runInfo.ecalPeriods))

#------------------------------------------------------------------------------
# Functions that are run for each file
#------------------------------------------------------------------------------

    def getCfgFileName(self, run, subrun):
        if self.runInfo.output == 'nd280':
            return self.runInfo.basePath+'/staging/nd280/cfg/nd280-{0}-{1}.cfg'.format(run,subrun)
        if self.runInfo.output == 'numc':
            return self.runInfo.basePath+'/staging/numc/cfg/numc-{0}-{1}.cfg'.format(run,subrun)
        if self.runInfo.output == 'nucp':
            return self.runInfo.basePath+'/staging/nucp/cfg/nucp-{0}-{1}.cfg'.format(run,subrun)


    def setInputFileList(self):
        self.inputFileList = []
        if self.runInfo.output == 'nd280':
            self.inputFileList.append(self.runInfo.PassThruDir+'/'+self.inputFile)
            if self.runInfo.startModule != 'nd280MC':
                # need to add numc file to inputFileList
                self.inputFileList.append(self.getNumcFile())
        else:
            # set inputFileList for nucp
            run = self.runInfo.RunNumber
            subrun = self.runInfo.SubRunNumber
            for a, b, files in os.walk(self.runInfo.PassThruDir):
                for aFile in files:
                    if '.root' in aFile:
                        if '%(run)03i-%(subrun)04i'%{'run':run%1000, 'subrun':subrun} in aFile:
                            self.inputFile = aFile
                            self.inputFileList.append(os.path.join(a,aFile))
            if len(self.inputFileList) == 0:
                print 'Could not find inputfiles', run, subrun
                self.runInfo.MissingConfigs.append((run,subrun))


    def GetShortPathFileList(self):
        string = ' '
        for file in self.inputFileList:
            sFile = file[1+file.rfind('/'):]
            string += ' '+sFile
        return string


    def GetFullPathFileList(self):
        string = ' '
        for file in self.inputFileList:
            string += ' '+file
        return string


    def getNumcFile(self):
        # If the startmodule is set as anything other than the default (nd280mc)
        # the inputFile will be different from the numc file.  Both are
        # required though, as the numc will be used by oaAnalysis for the
        # final passthrough and the inputFile will be passed to runND280 via
        # the .cfg file.    The critical point is to ensure that the numc file
        # is in the same directory as where nd280 was run, so this will be
        # copied over by the bash script.
        run = self.runInfo.RunNumber
        subrun = self.runInfo.SubRunNumber
        if self.NumcDirFiles.has_key(self.runInfo.Key):
            return self.NumcDirFiles[self.runInfo.Key]
        else:
            print 'Couldn\'t find numcfile!', run, subrun
            print 'Can\'t create script for', run, subrun


    def setInputFile(self):
        self.inputFile = None
        run = self.runInfo.RunNumber
        subrun = self.runInfo.SubRunNumber
        if self.PassThruDirFiles.has_key(self.runInfo.Key):
            self.inputFile = self.PassThruDirFiles[self.runInfo.Key].split('/')[-1]
        else:
            print 'Couldn\'t find inputfile!', run, subrun
            print 'Can\'t create script for', run, subrun


    def ResetRun(self):
        if '[configuration]' in self.cfgDict.keys():
            cfgList = self.cfgDict['[configuration]']
            for item in cfgList:
                if 'inputfile' in item:
                    cfgList.pop(cfgList.index(item))
        if '[filenaming]' in self.cfgDict.keys():
            fnList = self.cfgDict['[filenaming]']
            for item in fnList:
                if 'subrun' in item:
                    fnList.pop(fnList.index(item))
            for item in fnList:
                if 'run_number' in item:
                    fnList.pop(fnList.index(item))
        if '[neutrino]' in self.cfgDict.keys():
            nuList = self.cfgDict['[neutrino]']
            for item in nuList:
                if 'random_seed' in item:
                    nuList.pop(nuList.index(item))
        if '[electronics]' in self.cfgDict.keys():
            elList = self.cfgDict['[electronics]']
            for item in elList:
                if 'random_seed' in item:
                    elList.pop(elList.index(item))
        if '[nd280mc]' in self.cfgDict.keys():
            ndList = self.cfgDict['[nd280mc]']
            for item in ndList:
                if 'random_seed' in item:
                    ndList.pop(ndList.index(item))
        if '[cherry_picker]' in self.cfgDict.keys():
            cpList = self.cfgDict['[cherry_picker]']
            for item in cpList:
                if 'inputfile_list' in item:
                    cpList.pop(cpList.index(item))


    def setRunNumbersAndReset(self, run, subrun):
        self.ResetRun()
        run += self.runInfo.runprefix
        self.runInfo.RunNumber = run
        self.runInfo.SubRunNumber = subrun
        self.runInfo.Key = '{0}-{1:04}'.format(run,subrun)
        self.throwAgain()
        self.cfgDict['[filenaming]'].insert(0,('subrun', self.runInfo.SubRunNumber))
        self.cfgDict['[filenaming]'].insert(0,('run_number', self.runInfo.RunNumber))
        if self.runInfo.output == 'nd280':
            self.setInputFile()
            if self.inputFile:
                self.setInputFileList()
                self.cfgDict['[configuration]'].append(('inputfile', self.inputFile))
                self.cfgDict['[electronics]'].insert(0,('random_seed',random.randint(0,9e8)))
                self.cfgDict['[nd280mc]'].insert(0,('random_seed',random.randint(0,9e8)))
        if self.runInfo.output == 'numc':
            self.cfgDict['[neutrino]'].insert(0,('random_seed',random.randint(0,9e8)))
        if self.runInfo.output == 'nucp':
            self.setInputFileList()
            self.cfgDict['[nd280mc]'].insert(0,('random_seed',random.randint(0,9e8)))
            self.cfgDict['[cherry_picker]'].insert(0,('inputfile_list',self.GetShortPathFileList()))


    def setComment(self):
    #ver, vol, base, mcType):
        ver = int(self.runInfo.production)
        vol = self.runInfo.fluxMasterVolume.lower()
        base = self.runInfo.baseline
        beam = self.runInfo.beam
        fill = self.runInfo.fill
        comment = 'prod{0:03}{1}{2}{3}{4}'.format(ver, vol, base[0:4], base[5:7],fill)
        if vol == 'magnet':
            if beam=='run1':
                comment += 'a'
            elif beam =='run2':
                comment += 'b'
            elif beam =='run3' or beam =='run4':
                comment += 'c'
            elif beam =='run5':
                comment += 'd'
            else:
                print 'Baseline, fluxVolume, mcType configuration not compatible'
                print vol, base, mcType
                return 'NOTVALID'
        else:
            #fluxVolume is basket
            if 'beam' not in beam and 'nu' not in beam:
                comment += beam
        self.runInfo.Comment = comment
#------------------------------------------------------------------------------
# Construction of the config file
#------------------------------------------------------------------------------

    def getSubLines(self, sub):
        lines = []
        if sub in self.cfgDict.keys():
            list = self.cfgDict[sub]
            for key, val in list:
                lines.append(key+' = '+str(val))
        return lines


    def constructCfgFile(self):
        lines = []
        subs = ['[software]','[configuration]','[filenaming]','[analysis]']
        subs.extend(['[geometry]','[neutrino]','[nd280mc]', '[electronics]',
                                 '[cherry_picker]', '[dead_channels]'])

        for sub in subs:
            block = self.getSubLines(sub)
            if len(block) > 0:
                lines.append(sub)
                lines.extend(block)
                lines.append('')

        return lines


    def writeCfgFile(self, run, subrun):
        print self.getCfgFileName(run,subrun)
        outFile = open(self.getCfgFileName(run,subrun), 'w')
        self.setRunNumbersAndReset(run, subrun)
        lines = self.constructCfgFile()
        for line in lines:
            outFile.write(line+'\n')
        outFile.close()


    def parseRunSubrunArgs(self,arg):
        runs = []
        if '-' in arg:
            try:
                aSplit = [int(k) for k in arg.split('-')]
                runs = range(aSplit[0], aSplit[1]+1)
            except:
                print 'faulty argument:', arg
                sys.exit()
        else:
            runs = [int(k) for k in arg.split(',')]
        return runs


    def submitJob(self,run, subrun):
        queue_check(getpass.getuser())
        output = self.runInfo.output
        script = self.createScript(int(run), int(subrun))
        if not script:
            print 'Skipping job submission!'
            return
        print 'Created script', script
        handle = self.runInfo.getJobName(run,subrun)
        cmd = 'Qsub -e -l lnxfarm -N '
        subcmd = cmd + handle + ' -o '+output+'/logs/'+handle+'.out '+script
        print subcmd
        os.system(subcmd)


    def Usage(self):
        print 'Usage:'
        print 'to specify run/subrun number(s) (two arguments):'
        print '  ', 'an individual run'
        print '             ', sys.argv[0], 'runN subrunN'
        print '  ', 'a range:'
        print '             ', sys.argv[0], 'runM-runN subrunN'
        print '             ', sys.argv[0], 'runN subruM-subrunN'
        print '             ', sys.argv[0], 'runM-runN subruM-subrunN'
        print 'a pickled set (perhaps of missing runs) (one argument):'
        print '             ', sys.argv[0], 'Missing.pkl'
        print 'use the specified ranges as coded within this file (one argument):'
        print '             ', sys.argv[0], 'hardcoded'
        print 'display this usage'
        print '             ', sys.argv[0], 'help'


    def createScript(self, run, subrun):
        output = self.runInfo.output
        self.setRunNumbersAndReset(run,subrun)

        outFileName = output+'.'+self.runInfo.baseline+'.'+str(run)+'.'+str(subrun)+'.sh'
        outFileName = './'+output+'/scripts/'+outFileName
        outputFile = open(outFileName, 'w')
        HERE = get_cwd()
        STORE = (tdir.HADOOP_PROD+HERE.split('GNProd/', 1)[-1]
                 if self.runInfo.store == 'hadoop' else HERE)
        print 'output file storage area is:', STORE.replace('staging', '')
        out = []
        out.append('#!/bin/bash\n\n')
        out.append('export HERE='+HERE)
        out.append('export STORE='+STORE)
        out.append('cd ${TMPDIR}')
        out.append('echo "Copying file"')

        if (self.runInfo.output == 'nd280' and self.inputFile) or self.runInfo.output == 'nucp':
            out.append('cp -u '+self.GetFullPathFileList() + ' .')
        elif self.runInfo.output == 'numc':
            needFiles = self.runInfo.getNeededFiles()
            out.append('mkdir -p /sge-batch/GNfiles/genie_evtrate\n')
            for file in needFiles:
                src_file = os.path.join(self.runInfo.fluxDir, file)
                out.append('cp -u '+src_file+' /sge-batch/GNfiles/'+file+' & wait')
                # Need to check that the file has been completely copied
                # If we don't do this, and another job starts on the same
                # node and tries to open the partially copied flux file,
                # the ROOT file will be truncated and genie_mc will enter
                # an infinite loop of warnings, causing log-files to bloat.
                chek_str = ('while [ $(du -b /sge-batch/GNfiles/'+file+''
                            ' | cut -f 1) -lt $(du -b '+src_file+''
                            ' | cut -f 1) ]; do')
                out.append(chek_str)
                out.append('    sleep 5')
                out.append('done\n')
            for file in needFiles:
                out.append('ln -sf /sge-batch/GNfiles/'+file+' . &')
            out.append('wait\n')
        else:
            print 'Something went wrong with script creation!'
            return False

        out.append('\n')
        out.append('echo "Running ND280"')
        out.append('runND280 -t ${TMPDIR} -c '+self.getCfgFileName(run,subrun))
        out.append('\n')

        if self.runInfo.output == 'nd280':
            # store the big root files on hadoop
            to_store = ['cali', 'elmc', 'reco', 'anal', 'g4mc']
            # save log files on nfs
            to_save = ['logf']
        elif self.runInfo.output == 'numc':
            to_store = ['gnmc', 'numc']
            to_save = ['logf']
        elif self.runInfo.output == 'nucp':
            out.append('rm oa_gn_*nucp*geo.root')
            to_store = ['nucp']
            to_save = ['logf']

        for a_store in to_store:
            make_dirs_if_needed(STORE+'/../'+a_store)
            out.append('/usr/local/adm/bin/CUStageOut oa_gn*_'+a_store+'_* $STORE/../'+a_store)

        for a_save in to_save:
            make_dirs_if_needed(HERE+'/../'+a_save)
            out.append('/usr/local/adm/bin/CUStageOut oa_gn*_'+a_save+'_* $HERE/../'+a_save)

        for line in out:
            outputFile.write(line+'\n')
        outputFile.close()
        os.system('chmod 744 '+outFileName)
        return outFileName
