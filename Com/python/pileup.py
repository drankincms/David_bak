#!/usr/bin/env python
#########################################################################
#
# pileup.py
#
# pileup module
#
# Usage: ./analyze.py --mc-pileup <scenario>
#
#
# Authors:
#          Gena Kukartsev, July 2012
#
#
#########################################################################



import os
import glob
import re
import sys

from inspect import getmembers


legend = '[pileup]:'



def GetPuFiles(package):
    
    #get the release base location
    cmssw_release_base = os.environ['CMSSW_RELEASE_BASE']

    # path
    _path = cmssw_release_base+'/src/'+package+'/python'

    # get lists of py files
    _files = glob.glob(_path+'/*.py')

    return _files



def GetPuScenarios(files):

    _scenarios = []
    for path in files:
        _scenario = os.path.basename(path)
        _scenario = re.sub('_cff.py', '', _scenario)
        _scenario = re.sub('_cfi.py', '', _scenario)
        _scenario = re.sub('.py', '', _scenario)
        if _scenario.find('_') <= 0:
            continue
        # remove up to the first underscore
        #_scenario = _scenario[_scenario.find('_')+1:]
        _scenarios.append(_scenario)
        #print legend, _scenario

    return _scenarios



def GetMatchingScenarios(scenarios, pattern):
    #
    # find all PU scenarios that match pattern
    #
    _match = []
    for sc in scenarios:
        if sc.find(pattern)>=0:
            _match.append(sc)

    return _match



def RecursiveHasAttr(obj, attr):

    try:
        left, right = attr.split('.', 1)
    except:
        return hasattr(obj, attr)

    return RecursiveHasAttr(getattr(obj, left, None), right)



def RecursiveGetAttr(obj, attr, default=None):
    
    try:
        left, right = attr.split('.', 1)
    except:
        return getattr(obj, attr, default)
    
    return RecursiveGetAttr(getattr(obj, left, default), right, default)



def SavePuHist(files, sc, full_sim):
    #
    # saves the 'sc' scenario to a histogram
    #

    # possible names of the objects which contain PU info
    mod_names = []
    mod_names.append('famosPileUp.PileUpSimulator')
    mod_names.append('mix.input.nbPileupEvents')
    mod_names.append('mixGenPU.input.nbPileupEvents')

    _filename = 'pileup_'
    if full_sim:
        _filename += 'FULLSIM_'
    else:
        _filename += 'FASTSIM_'
    _filename = _filename + sc+'.root'


    _nfound = 0
    for file in files:
        if file.find(sc+'_cfi')>=0:
            _nfound += 1
            _module = re.sub('.py', '', os.path.basename(file))
            sys.path.append(os.path.dirname(file))

            exec('import '+_module+' as pu_mod')

            for mname in mod_names:
                if RecursiveHasAttr(pu_mod, mname):
                    _mod = RecursiveGetAttr(pu_mod, mname)

            _bin  = _mod.probFunctionVariable
            _prob = _mod.probValue

            _nbins = len(_bin)
            import ROOT
            _file = ROOT.TFile(_filename, 'recreate')
            _hist = ROOT.TH1F('MC_pileup', 'MC_pileup', _nbins, 0, _nbins)
            for bin in _bin:
                _hist.SetBinContent(bin+1, _prob[bin])

            print legend, 'saving to file', _filename
            _file.Write()

    if _nfound < 1:
        print legend, 'WARNING: PileUp information not found, cannot create a histogram!'

    return



    
def GetCvsHistogram(options):
    
    fastsim_cvs = 'FastSimulation/PileUpProducer'
    fullsim_cvs = 'SimGeneral/MixingModule'

    fastsim_files = GetPuFiles(fastsim_cvs)
    fullsim_files = GetPuFiles(fullsim_cvs)


    fastsim_scenarios = GetPuScenarios(fastsim_files)
    fullsim_scenarios = GetPuScenarios(fullsim_files)

    if options.mc_pileup != 'ALL':

        if options.fastsim_exact:
            fastsim_scenarios = [options.mc_pileup]
            fullsim_scenarios = []

        elif options.fullsim_exact:
            fullsim_scenarios = [options.mc_pileup]
            fastsim_scenarios = []

        else:
            fastsim_scenarios = GetMatchingScenarios(fastsim_scenarios,
                                                     options.mc_pileup)
            fullsim_scenarios = GetMatchingScenarios(fullsim_scenarios,
                                                     options.mc_pileup)


    if options.mc_pileup == 'FILES':
        for _file in fastsim_files:
            print legend, 'FASTSIM:', _file
        for _file in fullsim_files:
            print legend, 'FULLSIM:', _file


    if ( len(fastsim_scenarios) + len(fullsim_scenarios) ) < 1:
        print legend, 'no matching pileup scenarios found'
        return

    elif ( len(fastsim_scenarios) + len(fullsim_scenarios) ) > 1:
        print legend, 'more than one PU scenario found'
        print legend, 'possible FastSim pileup scenarios:'
        for sc in fastsim_scenarios:
            print legend, sc
        print legend, 'possible FullSim pileup scenarios:'
        for sc in fullsim_scenarios:
            print legend, sc

        return


    if len(fastsim_scenarios)==1:
        print legend, 'one FastSim pileup scenario found'
        SavePuHist(fastsim_files, fastsim_scenarios[0], full_sim = False)

    if len(fullsim_scenarios)==1:
        print legend, 'one FullSim pileup scenario found'
        SavePuHist(fullsim_files, fullsim_scenarios[0], full_sim = True)

    return
