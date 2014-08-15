import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes

import os

process = cms.Process("LJMetCom")

# Dilepton calculator options
process.load('LJMet.Com.DileptonCalc_cfi')

############################################################
#
# FWLite application options
#
process.ljmet = cms.PSet(
    isMc = cms.bool(True),
    runs = cms.vint32([])
)

process.DileptonCalc.isMc         = process.ljmet.isMc
process.DileptonCalc.dataType     = cms.string('MC')
process.DileptonCalc.genParticles = cms.InputTag("genParticles")

############################################################
#
# Event selector options
#
process.event_selector = cms.PSet(

    selection = cms.string('DileptonSelector'),

    # cuts
    #HLT
    trigger_cut              = cms.bool(True),
    dump_trigger             = cms.bool(True),

    #Can use same trigger paths for data and MC since MC is always one of the data versions
    trigger_path_ee          = cms.vstring('HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v15',
                               'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v16',
                               'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v17',
                               'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v18'),
    
    trigger_path_em          = cms.vstring('HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4', 'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5',
                               'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6', 'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7',
                               'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8', 'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9',
                               'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4', 'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5',
                               'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6', 'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7',
                               'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8', 'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9'),

    trigger_path_mm          = cms.vstring('HLT_Mu17_Mu8_v16', 'HLT_Mu17_Mu8_v17', 'HLT_Mu17_Mu8_v18', 'HLT_Mu17_Mu8_v19', 'HLT_Mu17_Mu8_v21'),    

    pv_cut                   = cms.bool(True),
    hbhe_cut                 = cms.bool(True),

    jet_cuts                 = cms.bool(False),
    jet_minpt                = cms.double(30.0),
    jet_maxeta               = cms.double(5),
    min_jet                  = cms.int32(0),
    max_jet                  = cms.int32(4000),

    muon_cuts                = cms.bool(True),
    min_muon                 = cms.int32(0),
    muon_minpt               = cms.double(10.0),
    muon_maxeta              = cms.double(2.5),
    max_muon                 = cms.int32(10),

    electron_cuts            = cms.bool(True),
    min_electron             = cms.int32(0),
    electron_minpt           = cms.double(10.0),
    electron_maxeta          = cms.double(2.5),
    max_electron             = cms.int32(10),

    min_lepton               = cms.int32(2),

    met_cuts                 = cms.bool(False),
    min_met                  = cms.double(0.0),
    btag_cuts                = cms.bool(False),
    btagger                  = cms.string('simpleSecondaryVertexHighEffBJetTags'),
    btag_min_discr           = cms.double(2.0),
    btag_1                   = cms.bool(True),
    btag_2                   = cms.bool(True),
    btag_3                   = cms.bool(False),

    trigger_collection       = cms.InputTag('TriggerResults::HLT'),
    pv_collection            = cms.InputTag('goodOfflinePrimaryVertices'),
    jet_collection           = cms.InputTag('goodPatJetsPFlow'),
    muon_collection          = cms.InputTag('selectedPatMuonsPFlow'),
    electron_collection      = cms.InputTag('selectedPatElectronsPFlow'),
    met_collection           = cms.InputTag('patMETsPFlow')

)

#######################################################
#
# Input files
#
input_module = 'LJMet.Com.Samples_2012.Dilepton.cffs.T53_950_cff'
process.load(input_module)

rel_base = os.environ['CMSSW_BASE']

# JSON
JsonFile = rel_base+'/src/LJMet/Com/data/json/Cert_190456-202016_8TeV_PromptReco_Collisions12_JSON_MuonPhys.txt'
myList   = LumiList.LumiList(filename=JsonFile).getCMSSWString().split(',')
if not process.ljmet.isMc:
    process.inputs.lumisToProcess.extend(myList)
        
        
        
#######################################################
#
# Output
#
process.outputs = cms.PSet (
    outputName = cms.string('T53T53_950'),
    treeName   = cms.string('ljmet'),
)

#######################################################
#
# Object selector options
#

# Primary vertex
process.load('PhysicsTools.SelectorUtils.pvSelector_cfi')
process.pvSelector.pvSrc   = cms.InputTag('goodOfflinePrimaryVertices')
process.pvSelector.minNdof = cms.double(4.0)
process.pvSelector.maxZ    = cms.double(24.0)
process.pvSelector.maxRho  = cms.double(2.0)

# Only requirements for muons: GlobalMuon and TrackerMuon
process.load('PhysicsTools.SelectorUtils.pfMuonSelector_cfi') 
process.pfMuonSelector.cutsToIgnore = cms.vstring(   'Chi2',
                                                     'D0',
                                                     'NHits',
                                                     'NValMuHits',
                                                     'PFIso',
                                                     'nPixelHits',
                                                     'nMatchedStations')

# Most electron requirements in place, but not isolation or d0
process.load('LJMet.Com.cutbasedIDSelector_cfi')
process.cutbasedIDSelector.cutsToIgnore = cms.vstring('reliso_EB',
                                                      'reliso_EE',
                                                      'd0_EB',
                                                      'd0_EE')

# jets
process.load('PhysicsTools.SelectorUtils.pfJetIDSelector_cfi') 
process.pfJetIDSelector.version = cms.string('FIRSTDATA')
process.pfJetIDSelector.quality = cms.string('LOOSE')
