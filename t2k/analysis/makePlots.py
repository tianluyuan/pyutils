"""
Modified from truthInfoPlot.py to plot stacked histos with true
reaction code info.
T.Y. 6/24/2012  Modifiled T.T. 7/2/2012, 9/5/2012
Moved to Python/tianlu/analysis module 2/23/2014 for global usage

USAGE (deprecated):
python makePlots.py pathTo/dataFile.root pathTo/mcFile.root
variableToHist binLow binHigh binStep setLogybool histOption fluxRatio
cutOption optTitle
"""

import os
import glob
import numpy as np
import ROOT as rt
from array import array
from tianlu import directories
from tianlu.analysis.P0DCutUtils import (getPOTsFromChain,
                                         goodColorIndices,
                                         prompt)
from tianlu.rootutils import get_legend, glob_to_chain
from operator import itemgetter


__author__ = 'Tianlu Yuan'
__email__ = 'tianlu.yuan [at] colorado.edu'


class makeComparisonPlots:
    def __init__(self, filePath_data, filePath_mc):
        # Data and MC filepaths
        self._filePath_data = filePath_data
        self._filePath_mc = filePath_mc

        self.dataFilesName = os.path.relpath(
            os.path.dirname(self._filePath_data),
            directories.RUNANALYSES_STORE).replace('/', '_')
        self.mcFilesName = os.path.relpath(
            os.path.dirname(self._filePath_mc),
            directories.RUNANALYSES_STORE).replace('/', '_')

        self.ch1_data = rt.TChain('t1', 'flatTreeChain')
        self.ch1_mc = rt.TChain('t1', 'flatTreeChain')
        self.ch2_data = rt.TChain('t2', 'counterTreeChain_data')
        self.ch2_mc = rt.TChain('t2', 'counterTreeChain_mc')

        glob_to_chain(self.ch1_data, self._filePath_data)
        glob_to_chain(self.ch1_mc, self._filePath_mc)

        glob_to_chain(self.ch2_data, self._filePath_data)
        glob_to_chain(self.ch2_mc, self._filePath_mc)

        self.pot_data = getPOTsFromChain(self.ch2_data)[0]
        # To normalize combined beam/sand mc correctly we need to keep
        # track of their separate pots.
        self.pot_mc, self.pot_mc_sand = getPOTsFromChain(self.ch2_mc)
        print ('self.pot_data:', self.pot_data,
               ', self.pot_mc:', self.pot_mc,
               ', self.pot_mc_sand:', self.pot_mc_sand)

        # Normalize by POT
        # dont wanna divide by zero
        try:
            self.data_pot_inv = str(1/self.pot_data)
        except ZeroDivisionError:
            print 'pot_data is 0!'
            self.data_pot_inv = str(0)
        try:
            self.mc_pot_inv = str(1/self.pot_mc)
        except ZeroDivisionError:
            print 'pot_mc is 0!'
            self.mc_pot_inv = str(0)
        try:
            self.mc_sand_pot_inv = str(1/self.pot_mc_sand)
        except ZeroDivisionError:
            self.mc_sand_pot_inv = str(0)

        # Set up canvases
        self._myCanvas = rt.TCanvas('myCanvas', 'My Canvas')
        self._c1 = rt.TCanvas('c1')


    def setupVarsFromInput(self, outdir, **kwarg):

        production = self._filePath_data.split('prod')[1][0:2]

        # Info on variable to be histogrammed
        self.variableName = kwarg['varname']
        self.variableToPlot = kwarg['varrootname']

        # Output plot file directory
        self.outfiledir = outdir + self.variableName + '/'

        # make directory if it doesn't already exist
        if not os.path.exists(self.outfiledir):
            os.mkdir(self.outfiledir)

    #     # Charge selection
    #     charge = sys.argv[5]
    #     if int(charge) == -1:
    #         chargeStr = 'negCharge'
    #     elif int(charge) == 1:
    #         chargeStr = 'posCharge'
    #     elif int(charge) == 0:
    #         chargeStr = 'allCharge'
    #     else:
    #         print 'Bad charge inputted!'
    #         sys.exit(0)

        # binning requirements
        if kwarg.has_key('varbins'):
            self.binLowX = kwarg['varbins']
        else:
            binLow = kwarg['binlow']
            binHigh = kwarg['binhigh']
            binStep = kwarg['binstep']
            self.binLowX = np.arange(binLow, binHigh, binStep)

        # Set Logy axis 0 or 1
        self.setLogybool = bool(kwarg['logbool'])

        # normalization options
        self.normStr = 'norm2POT'
        self.yaxisTit = 'Count/POT'

        # Set norm to unity bool for Data/MC comparison
        self.normToUnityBool = bool(int(kwarg['norm']))
        if self.normToUnityBool:
            self.normStr = 'norm2Unity'
            self.yaxisTit = 'Density'

        # option to reweight MC by tuned flux 11bv3.2
        self.fluxStr = kwarg['tunedFlx']
        if self.fluxStr != 'FluxRatio':
            # if we dont want to used tuned flux, treat as if
            # we're reweighting events by unity
            self.fluxStr = '1'
            self.outfileFluxTitle = ''
        else:
            self.outfileFluxTitle = 'tuned11b32_'

        # any additional draw/cut options
        try:
            self.optionalCuts = kwarg['optCut']
            self.optCutName = kwarg['optCutName']
            self.stackTitle = self.normStr+', '+self.optCutName
            print 'Plotting with optional cuts: ' +self.optionalCuts
        except KeyError:
            self.optionalCuts = '1'
            self.optCutName = ''
            self.stackTitle = self.normStr+', noopt'

        # output options
        # root hists file
        self.rootHistosFile = rt.TFile(self.outfiledir+'../'+kwarg['rootHistosFile'], 'UPDATE')
        self.printCanvas = kwarg['printCanvas']


    def plot(self, outdir, **kwarg):
        """ Template for plotting
        """
        self.setupVarsFromInput(outdir, **kwarg)

        # These will be histos of the energy distribution
        # MC hists
        hist_ccqes = rt.TH1F('hist_ccqes', '#nu_{#mu} CCQE',
                          len(self.binLowX)-1, array('f', self.binLowX))
        hist_ccqes.Sumw2()
        hist_ccqes.GetXaxis().SetTitle(self.variableName)
        hist_ccqes.GetYaxis().SetTitle(self.yaxisTit)

        hist_ccres = rt.TH1F('hist_ccres', '#nu_{#mu} CCRes',
                          len(self.binLowX)-1, array('f', self.binLowX))
        hist_ccres.Sumw2()
        hist_ccres.GetXaxis().SetTitle(self.variableName)
        hist_ccres.GetYaxis().SetTitle(self.yaxisTit)

        hist_ccoth = rt.TH1F('hist_ccoth', 'CCOther',
                          len(self.binLowX)-1, array('f', self.binLowX))
        hist_ccoth.Sumw2()
        hist_ccoth.GetXaxis().SetTitle(self.variableName)
        hist_ccoth.GetYaxis().SetTitle(self.yaxisTit)

        hist_ncall = rt.TH1F('hist_ncall', 'NC',
                          len(self.binLowX)-1, array('f', self.binLowX))
        hist_ncall.Sumw2()
        hist_ncall.GetXaxis().SetTitle(self.variableName)
        hist_ncall.GetYaxis().SetTitle(self.yaxisTit)

        hist_outfv = rt.TH1F('hist_outfv', 'Out P0D FV',
                          len(self.binLowX)-1, array('f', self.binLowX))
        hist_outfv.Sumw2()
        hist_outfv.GetXaxis().SetTitle(self.variableName)
        hist_outfv.GetYaxis().SetTitle(self.yaxisTit)

        hist_notru = rt.TH1F('hist_notru', 'No Truth',
                          len(self.binLowX)-1, array('f', self.binLowX))
        hist_notru.Sumw2()
        hist_notru.GetXaxis().SetTitle(self.variableName)
        hist_notru.GetYaxis().SetTitle(self.yaxisTit)

        hist_other = rt.TH1F('hist_other', 'Other',
                          len(self.binLowX)-1, array('f', self.binLowX))
        hist_other.Sumw2()
        hist_other.GetXaxis().SetTitle(self.variableName)
        hist_other.GetYaxis().SetTitle(self.yaxisTit)

        # Data hist
        hist_data = rt.TH1F('hist_data', 'Data',
                         len(self.binLowX)-1, array('f', self.binLowX))
        hist_data.Sumw2()
        hist_data.GetXaxis().SetTitle(self.variableName)
        hist_data.GetYaxis().SetTitle(self.yaxisTit)

        self._c1.cd()
        mc_base_draw_selection = (self.fluxStr+'*(!IsSandMC*'+
                                  self.mc_pot_inv+'+IsSandMC*'+
                                  self.mc_sand_pot_inv+')'+
                                  '*(%(stack)s'+
                                  self.optionalCuts+')')
        # Selecting events based on optional cuts
        self.ch1_mc.Draw(self.variableToPlot+'>>hist_ccqes',
                         mc_base_draw_selection %
                         dict(stack=('(trueGlobalMatch&&vtxGlobalMatch)&&'
                                     'TrueVtxInP0DFiducial&&'
                                     'ReactionInt % 100 == 11&&'
                                     'TruePDG == 13&&')))
        self.ch1_mc.Draw(self.variableToPlot+'>>hist_ccres',
                         mc_base_draw_selection %
                         dict(stack=('(trueGlobalMatch&&vtxGlobalMatch)&&'
                                     'TrueVtxInP0DFiducial&&'
                                     'ReactionInt % 100 == 12&&'
                                     'TruePDG == 13&&')))
        self.ch1_mc.Draw(self.variableToPlot+'>>hist_ccoth',
                         mc_base_draw_selection %
                         dict(stack=('(trueGlobalMatch&&vtxGlobalMatch)&&'
                                     '(((ReactionInt % 100 == 11||'
                                     'ReactionInt % 100 == 12)&&'
                                     'TruePDG != 13)||'
                                     'ReactionInt % 100 == 13 ||'
                                     'ReactionInt % 100 == 14)&&'
                                     'TrueVtxInP0DFiducial&&')))
        self.ch1_mc.Draw(self.variableToPlot+'>>hist_ncall',
                         mc_base_draw_selection %
                         dict(stack=('(trueGlobalMatch&&vtxGlobalMatch)&&'
                                     'TrueVtxInP0DFiducial&&'
                                     '(ReactionInt % 100 > 20)&&')))
        self.ch1_mc.Draw(self.variableToPlot+'>>hist_outfv',
                         mc_base_draw_selection %
                         dict(stack=('(trueGlobalMatch&&vtxGlobalMatch)&&'
                                     '!TrueVtxInP0DFiducial&&')))
        self.ch1_mc.Draw(self.variableToPlot+'>>hist_notru',
                         mc_base_draw_selection %
                         dict(stack=('!(trueGlobalMatch&&vtxGlobalMatch)&&')))
        self.ch1_mc.Draw(self.variableToPlot+'>>hist_other',
                         mc_base_draw_selection % dict(stack=''))
        self.ch1_data.Draw(self.variableToPlot+'>>hist_data',
                           self.data_pot_inv+'*('+self.optionalCuts+')')

        # For some reason ROOT or pyROOT forces you to set fill color after
        # hist has been filled
        hist_ccqes.SetFillColor(2)
        hist_ccqes.SetLineColor(2)
        hist_ccres.SetFillColor(5)
        hist_ccres.SetLineColor(5)
        hist_ccoth.SetFillColor(3)
        hist_ccoth.SetLineColor(3)
        hist_ncall.SetFillColor(7)
        hist_ncall.SetLineColor(7)
        hist_outfv.SetFillColor(6)
        hist_outfv.SetLineColor(6)
        hist_notru.SetFillColor(4)
        hist_notru.SetLineColor(4)
        hist_other.SetFillColor(12)
        hist_other.SetLineColor(12)

        # Legends
        # legPlacement = (0.35, 0.45, 0.68, 0.88)
        legTextSize = 0.04
        leg = get_legend()
        leg.SetTextSize(legTextSize)
        leg.AddEntry(hist_ccqes, hist_ccqes.GetTitle(), 'f')
        leg.AddEntry(hist_ccres, hist_ccres.GetTitle(), 'f')
        leg.AddEntry(hist_ccoth, hist_ccoth.GetTitle(), 'f')
        leg.AddEntry(hist_ncall, hist_ncall.GetTitle(), 'f')
        leg.AddEntry(hist_outfv, hist_outfv.GetTitle(), 'f')
        leg.AddEntry(hist_notru, hist_notru.GetTitle(), 'f')
        leg.AddEntry(hist_other, hist_other.GetTitle(), 'f')
        rt.SetOwnership(leg, 0)

        # Stack histos
        stack = rt.THStack('stack', self.stackTitle)
        stack.Add(hist_notru)
        stack.Add(hist_outfv)
        stack.Add(hist_ncall)
        stack.Add(hist_ccoth)
        stack.Add(hist_ccres)
        stack.Add(hist_ccqes)

        self._myCanvas.cd()
        if self.setLogybool:
            self._myCanvas.SetLogy()
        else:
            self._myCanvas.SetLogy(0)

        # Not sure why but need 'hist' option to draw filled histos
        # The following four histos are saved in histos.root file
        # in the appropriate plots directory.  Their titles
        # are set so that they make sense in the root file.
        hist_other.SetTitle('mc '+self.stackTitle.split(', ')[1])
        hist_data.SetTitle('data '+self.stackTitle.split(', ')[1])

        if self.normToUnityBool:
            hist_other.SetFillColor(0)
            hist_other.SetLineColor(2)
            hist_other.DrawNormalized('histE')
            hist_data.DrawNormalized('sameE')
        else:
            # Find max between hist_other and hist_data and sets as max for both
            mc_max = hist_other.GetMaximum()
            data_max = hist_data.GetMaximum()
            if mc_max > data_max:
                hist_data.SetMaximum(mc_max * 1.1)
            else:
                hist_other.SetMaximum(data_max * 1.1)

            hist_other.Draw('hist')
            stack.Draw('histsame')
            hist_data.Draw('sameE')
            # Must call draw first before setting rt.THStack axis title 
            stack.GetXaxis().SetTitle(self.variableName)
            leg.Draw()

        # save mc and data histos in root file
        #
        # overwrite if existing key exists.  only write if non-empty
        hist_name_base = '{var}_{flux_title}{data_sample}_'\
                         '{normalization}{option_cut}'
        template = hist_name_base.format(var = self.variableName,
                                         flux_title = '{flux_title}',
                                         data_sample = '{data_sample}',
                                         normalization = self.normStr,
                                         option_cut = self.optCutName,)

        if hist_other.GetEntries() > 0:
            hist_other_name = template.format(flux_title = self.outfileFluxTitle,
                                              data_sample = self.mcFilesName,)
            hist_other.Write(hist_other_name, rt.TObject.kOverwrite)
            stack.Write(hist_other_name+'_thstack', rt.TObject.kOverwrite)
        if hist_data.GetEntries() > 0:
            hist_data_name = template.format(flux_title = '',
                                             data_sample = self.dataFilesName,)
            hist_data.Write(hist_data_name, rt.TObject.kOverwrite)

        canvas_name = template.format(flux_title = self.outfileFluxTitle,
                                      data_sample = self.dataFilesName+'_'\
                                      +self.mcFilesName,)

        self._myCanvas.Write(canvas_name+'_tcanvas', rt.TObject.kOverwrite)

        if self.printCanvas:
            canvas_filename_template = '{outfiledir}{flux_title}{data}'\
                                       '_{mc}_{normalization}_{option}'
            canvas_filename = canvas_filename_template.format(
                outfiledir = self.outfiledir,
                flux_title = self.outfileFluxTitle,
                data = self.dataFilesName,
                mc = self.mcFilesName,
                normalization = self.normStr,
                option = self.optCutName)

            self._myCanvas.Print(canvas_filename+'.pdf')
            self._myCanvas.Print(canvas_filename+'.eps')
            self._myCanvas.Print(canvas_filename+'.C')

        self.rootHistosFile.Close()


    def compareMCFluxTuning(self, **kwarg):
        self.setupVarsFromInput(**kwarg)

        # MC hist, no flux reweighting
        # used for MC flux rewieghting comparison
        hist_mcnotune = rt.TH1F('hist_mcnotune', 'Other', len(self.binLowX)-1, array('f', self.binLowX))
        hist_mcnotune.Sumw2()
        hist_mcnotune.GetXaxis().SetTitle(self.variableName)
        hist_mcnotune.GetYaxis().SetTitle(self.yaxisTit)

        hist_mctune = rt.TH1F('hist_mctune', 'Other', len(self.binLowX)-1, array('f', self.binLowX))
        hist_mctune.Sumw2()
        hist_mctune.GetXaxis().SetTitle(self.variableName)
        hist_mctune.GetYaxis().SetTitle(self.yaxisTit)

        # For comparing MC to itself, i.e. tuned vs untuned flux
        self.ch1_mc.Draw(self.variableToPlot+'>>hist_mcnotune', '(IsMC &&'+self.optionalCuts+')')
        self.ch1_mc.Draw(self.variableToPlot+'>>hist_mctune', 'FluxRatio*(IsMC &&'+self.optionalCuts+')')
        
        hist_mcnotune.SetTitle('mcNoFlxTune '+self.stackTitle.split(', ')[1])
        hist_mctune.SetTitle('mcFlxTune '+self.stackTitle.split(', ')[1])

        hist_mctune.SetLineColor(2)
        hist_mctune.Draw('histE')
        hist_mcnotune.Draw('samehistE')

        # save to root hists
        # hist_mctune.Write(self.variableName+'_mcFlxTune_'+self.mcFilesName+'_'+self.optCutName)
        # hist_mcnotune.Write(self.variableName+'_mcNoFlxTune_'+self.mcFilesName+'_'+self.optCutName)

        if self.printCanvas:
            self._myCanvas.Print(self.outfiledir+'mcFlxTuneComp_'+self.mcFilesName+'_'+self.optCutName+'.pdf')
            self._myCanvas.Print(self.outfiledir+'mcFlxTuneComp_'+self.mcFilesName+'_'+self.optCutName+'.eps')
            self._myCanvas.Print(self.outfiledir+'mcFlxTuneComp_'+self.mcFilesName+'_'+self.optCutName+'.C')

        self.rootHistosFile.Close()

# make plots using the t3 tree
# T.Y. 6/26/13
class makePlotsFromT3:
    def __init__(self):
        _prod_list = ['prod5F']
        _sel_list = ['GetSandEvents','GetEvents','GetCosmicsEvents']        
 
        print 'Initial setup...'
        self._production = prompt('Please select a production from '+str(_prod_list)+': ',
                                  _prod_list)
        self._selection = prompt('Please select a selection from '+str(_sel_list)+': ',
                                 _sel_list)
        self._runlist = set(['rdpRUN3c_','mcpneutRUN3_', 'mcpneutRUN3sand_'])
        print 'Runlist is set as:',self._runlist

        self._full_run_set = set()
        self.getUniqueRuns()

        # chain_dict: {'some_RUN':(corresponding_Chain, color_index)}
        self._chain_dict = {}
        self.createChains()


    def getUniqueRuns(self):
        glob_str= 'outputs/%(selection)s_%(production)s_*.root'
        out_file_list = glob.glob(glob_str % dict(selection=self._selection, production=self._production))
        for a_file in out_file_list:
            a_run = a_file.split('_')[-2]+'_'
            self._full_run_set.add(a_run)


    def createChains(self):
        base_output_string = 'outputs/%(selection)s_%(production)s_%(run)s*.root'
        for idx, a_RUN in enumerate(self._runlist):
            self._chain_dict[a_RUN] = (rt.TChain('t3'), goodColorIndices[idx])

            filepath_to_chain = base_output_string % dict(selection=self._selection,
                                                          production=self._production,
                                                          run=a_RUN)
            print 'Adding',filepath_to_chain,'to',a_RUN,'chain'
            glob_to_chain(self._chain_dict[a_RUN][0], filepath_to_chain)


    def refresh(self):
        self.createChains()


    def appendToRunList(self):
        run_to_append = prompt('What run to append from '+str(self._full_run_set)+'? ',
                               self._full_run_set)

        if not run_to_append: return

        if run_to_append not in self._runlist: self._runlist.add(run_to_append)
        # refresh the chains
        self.refresh()


    def plotNPIDsPerBunch(self):
        variable = 'PIDsPerBunchArray'
        optional_cut = raw_input('Any optional cuts to apply? ')

        # hists_dict: {'hist_key':hist_object, ...}
        hists_dict = {}
        # ordered_hists: [(hist_max, hist_object), ...]
        ordered_hists = []
        leg = get_legend()
        leg.SetTextSize(0.04)
        
        for key, chain_tuple in self._chain_dict.iteritems():
            chain = chain_tuple[0]
            color = chain_tuple[1]

            hists_dict[key] = rt.TH1F('hist_'+key, 'PIDs Per Bunch', 20, 0, 20) 
            hists_dict[key].Sumw2()
            hists_dict[key].GetXaxis().SetTitle(variable)
            hists_dict[key].GetYaxis().SetTitle('Density')
            chain.Draw(variable+'>>hist_'+key, optional_cut)

            hists_dict[key].SetLineColor(color)

            hist_max = hists_dict[key].GetMaximum()
            hist_entries = hists_dict[key].GetEntries()
            ordered_hists.append((hist_max, hist_max/hist_entries, hists_dict[key]))

            leg.AddEntry(hists_dict[key], key, 'l')
        rt.SetOwnership(leg, 0)

        myCanvas = rt.TCanvas('myCanvas')
        myCanvas.cd()
        
        clear_canvas = True

        # sort hists by their maximum/entries since we are drawing normalized.  
        # Want the highest hist to be drawn first
        ordered_hists.sort(key=itemgetter(1), reverse=True)
        for hist_tuple in ordered_hists:
            hist_object = hist_tuple[2]
            draw_option = 'E'
            if 'mcp' in hist_object.GetName():
                draw_option+='hist'
            if not clear_canvas:
                draw_option+='same'

            hist_object.DrawNormalized(draw_option)
            clear_canvas = False

        leg.Draw()
        myCanvas.WaitPrimitive()

