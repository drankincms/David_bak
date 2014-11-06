import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes

process = cms.Process("LJMetCom")

############################################################
#
# FWLite application options
process.load('LJMet.Com.ljmet_cfi')
process.ljmet.isMc = cms.bool(False)

# common calculator options
process.load('LJMet.Com.commonCalc_cfi')
process.CommonCalc.dummy_parameter = cms.string('Dummy parameter value')

# Stop calculator options
process.load('LJMet.Com.stopCalc_cfi')

# Wprime calculator options
process.load('LJMet.Com.wprimeCalc_cfi')
process.WprimeCalc.isWJets = cms.bool(False)
process.WprimeCalc.isTB = cms.bool(False)
process.WprimeCalc.isTT = cms.bool(False)

# LjetsTopoCalc options
process.load('LJMet.Com.ljetsTopoCalcNew_cfi')
process.LjetsTopoCalcNew.useBestTop = cms.bool(True)


############################################################
#
# Event selector options
#
process.event_selector = cms.PSet(

    selection = cms.string('WprimeSelector'),

    # cuts
    debug  = cms.bool(False),
    isMc  = cms.bool(False),
    
    trigger_cut  = cms.bool(True),
    dump_trigger = cms.bool(False),
    
    mctrigger_path_el = cms.string('HLT_Ele27_WP80_v10'), 
    mctrigger_path_mu = cms.string('HLT_IsoMu24_eta2p1_v13'), 
    trigger_path_el = cms.vstring('HLT_Ele27_WP80_v8','HLT_Ele27_WP80_v9','HLT_Ele27_WP80_v10','HLT_Ele27_WP80_v11'), 
    trigger_path_mu = cms.vstring('HLT_IsoMu24_eta2p1_v11','HLT_IsoMu24_eta2p1_v12','HLT_IsoMu24_eta2p1_v13','HLT_IsoMu24_eta2p1_v14','HLT_IsoMu24_eta2p1_v15'),
 
    pv_cut         = cms.bool(True),
    hbhe_cut       = cms.bool(True),
    doLaserCalFilt = cms.bool(False),

    jet_cuts                 = cms.bool(True),
    jet_minpt                = cms.double(30.0),
    jet_maxeta               = cms.double(2.4),
    min_jet                  = cms.int32(2),
    max_jet                  = cms.int32(4000),
    leading_jet_pt           = cms.double(80.0),

    muon_cuts                = cms.bool(True),
    tight_muon_minpt         = cms.double(26.0), # 26 GeV
    tight_muon_maxeta        = cms.double(2.1),
    loose_muon_minpt         = cms.double(10.0),
    loose_muon_maxeta        = cms.double(2.4),
    min_tight_muon           = cms.int32(1),

    electron_cuts                = cms.bool(True),
    tight_electron_minpt         = cms.double(30.0), # 30 GeV
    tight_electron_maxeta        = cms.double(2.5),
    loose_electron_minpt         = cms.double(20.0),
    loose_electron_maxeta        = cms.double(2.5),
    min_tight_electron           = cms.int32(0),
    
    min_tight_lepton         = cms.int32(1),
    max_tight_lepton         = cms.int32(1),
    trigger_consistent       = cms.bool(True),
    second_lepton_veto       = cms.bool(True),
    
    met_cuts                 = cms.bool(True),
    min_met                  = cms.double(20.0),
    btag_cuts                = cms.bool(False),
    
    btagger                  = cms.string('combinedSecondaryVertexBJetTags'),
    btag_min_discr           = cms.double(0.679),
    btag_1                   = cms.bool(False),
    btag_2                   = cms.bool(False),
    btag_3                   = cms.bool(False),

    BTagUncertUp             = cms.bool(False),
    BTagUncertDown           = cms.bool(False),
    JECup                    = cms.bool(False),
    JECdown                  = cms.bool(False),
    JERup                    = cms.bool(False),
    JERdown                  = cms.bool(False),
    JEC_txtfile = cms.string('../cond/Fall12_V7_DATA_UncertaintySources_AK5PFchs.txt'),
    do53xJEC                 = cms.bool(True),
    trigger_collection       = cms.InputTag('TriggerResults::HLT'),
    pv_collection            = cms.InputTag('goodOfflinePrimaryVertices'),
    jet_collection           = cms.InputTag('goodPatJetsPFlow'),
    muon_collection          = cms.InputTag('selectedPatMuonsPFlow'),
    electron_collection      = cms.InputTag('selectedPatElectronsPFlow'),
    met_collection           = cms.InputTag('patMETsPFlow'),
    type1corrmet_collection  = cms.InputTag('pfType1CorrectedMet'),

    )


#######################################################
#
# Input files
#
process.inputs = cms.PSet (
   nEvents    = cms.int32(-1),
    skipEvents = cms.int32(0),
    lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange()),
    fileNames  = cms.vstring(
                           'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingleMu/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/SingleMu/USER/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/0000/F6C66144-482F-E211-835C-002618943836.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingleMu/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/SingleMu/USER/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/0000/F695947A-2E2F-E211-9997-003048FFD796.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingleMu/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/SingleMu/USER/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/0000/F68F0429-482F-E211-A530-002618FDA204.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingleMu/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/SingleMu/USER/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/0000/F68C64A9-472F-E211-93C2-003048678FA6.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingleMu/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/SingleMu/USER/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/0000/F688BB4B-482F-E211-B605-002618943877.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingleMu/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/SingleMu/USER/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/0000/F67CC017-472F-E211-AD65-001A92810AAE.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingleMu/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/SingleMu/USER/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/0000/F626B663-472F-E211-8BB6-003048FFD732.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingleMu/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/SingleMu/USER/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/0000/F604A1CC-2E2F-E211-99FD-002618943966.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingleMu/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/SingleMu/USER/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/0000/F4E0D0F5-2D2F-E211-8C47-003048678B12.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingleMu/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/SingleMu/USER/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/0000/F4D6BA38-472F-E211-B187-0026189438B8.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingleMu/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/SingleMu/USER/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/0000/F4BC222B-482F-E211-A410-002618943880.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingleMu/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/SingleMu/USER/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/0000/F4A49621-462F-E211-A312-001A92811744.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingleMu/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/SingleMu/USER/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/0000/F4629BA8-472F-E211-800E-00261894397E.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingleMu/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/SingleMu/USER/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/0000/F451884E-482F-E211-9609-0026189438AB.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingleMu/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/SingleMu/USER/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/0000/F2FF4514-2D2F-E211-81AD-002618943974.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingleMu/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/SingleMu/USER/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/0000/F2EC3702-472F-E211-AE71-003048679294.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingleMu/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/SingleMu/USER/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/0000/F2E68902-462F-E211-B8CD-003048678DA2.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingleMu/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/SingleMu/USER/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/0000/F2E2DB1D-462F-E211-8C8B-003048678BB8.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingleMu/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/SingleMu/USER/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/0000/F2D28488-2E2F-E211-BA2B-003048678B12.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingleMu/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/SingleMu/USER/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/0000/F2B03832-482F-E211-B475-002618943934.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingleMu/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/SingleMu/USER/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/0000/F2A25732-462F-E211-AC51-0026189438D6.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingleMu/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/SingleMu/USER/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/0000/F296EE89-452F-E211-B2F4-003048FFD732.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingleMu/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/SingleMu/USER/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/0000/F2606D12-2E2F-E211-950D-002618943966.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingleMu/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/SingleMu/USER/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/0000/F2450801-472F-E211-B824-002354EF3BDC.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingleMu/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/SingleMu/USER/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/0000/F209C6B5-472F-E211-806A-001A92971BC8.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingleMu/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/SingleMu/USER/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/0000/F0FBCA4C-2F2F-E211-AAD1-003048FFD71A.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingleMu/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/SingleMu/USER/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/0000/F0EAEC8E-452F-E211-9E52-001A92971B20.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingleMu/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/SingleMu/USER/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/0000/F0E056EB-472F-E211-B3CF-0026189438DA.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingleMu/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/SingleMu/USER/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/0000/F08CB131-462F-E211-BF0D-003048679266.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingleMu/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/SingleMu/USER/StoreResults-Run2012D-PromptReco-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/0000/F024CEAD-472F-E211-8DA0-003048678E80.root',

     )
)


# JSON
JsonFile = '../data/json/Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON_203894-208686.txt'
myList   = LumiList.LumiList(filename=JsonFile).getCMSSWString().split(',')
if (not process.ljmet.isMc == cms.bool(True)):
    process.inputs.lumisToProcess.extend(myList)
        
#######################################################
#
# Output
#
process.outputs = cms.PSet (
    outputName = cms.string('SingleMu_Prompt2012D1_2'),
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

# Tight muon
process.load('PhysicsTools.SelectorUtils.pfMuonSelector_cfi')
process.pfMuonSelector.version = cms.string('TOPPAG12_LJETS')

# Loose muon                                                       
process.looseMuonSelector = process.pfMuonSelector.clone()
process.looseMuonSelector.version = cms.string('TOPPAG12_LJETS_VETO')

# electron
process.load('LJMet.Com.cutbasedIDSelector_cfi')
process.cutbasedIDSelector.version = cms.string('TIGHT')
process.cutbasedIDSelector.cutsToIgnore     = cms.vstring()

# loose electron
process.looseElectronSelector = process.cutbasedIDSelector.clone()
process.looseElectronSelector.version = cms.string('VETO')
process.looseElectronSelector.cutsToIgnore     = cms.vstring()

# jets
process.load('PhysicsTools.SelectorUtils.pfJetIDSelector_cfi') 
process.pfJetIDSelector.version = cms.string('FIRSTDATA')
process.pfJetIDSelector.quality = cms.string('LOOSE')

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("histo.root"),
                                   closeFileFast = cms.untracked.bool(True)
                                   )

  
