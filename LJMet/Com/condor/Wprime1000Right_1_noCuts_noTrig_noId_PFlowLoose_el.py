import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes

process = cms.Process("LJMetCom")

############################################################
#
# FWLite application options
process.load('LJMet.Com.ljmet_cfi')
process.ljmet.isMc = cms.bool(True)

# common calculator options
process.load('LJMet.Com.commonCalc_cfi')
process.CommonCalc.dummy_parameter = cms.string('Dummy parameter value')

# Stop calculator options
process.load('LJMet.Com.stopCalc_cfi')

# Wprime calculator options
process.load('LJMet.Com.wprimeCalc_cfi')
process.WprimeCalc.isWJets = cms.bool(False)
process.WprimeCalc.isTB = cms.bool(True)
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
    isMc  = cms.bool(True),
    
    trigger_cut  = cms.bool(False),
    dump_trigger = cms.bool(False),
    
    mctrigger_path_el = cms.string('HLT_Ele27_WP80_v10'), 
    #mctrigger_path_mu = cms.string('HLT_IsoMu24_eta2p1_v13'),
    mctrigger_path_mu = cms.string('HLT_Ele27_WP80_v10'),
    trigger_path_el = cms.vstring('HLT_Ele27_WP80_v8','HLT_Ele27_WP80_v9','HLT_Ele27_WP80_v10','HLT_Ele27_WP80_v11'), 
    trigger_path_mu = cms.vstring('HLT_IsoMu24_eta2p1_v11','HLT_IsoMu24_eta2p1_v12','HLT_IsoMu24_eta2p1_v13','HLT_IsoMu24_eta2p1_v14','HLT_IsoMu24_eta2p1_v15'),
 
    pv_cut         = cms.bool(True),
    hbhe_cut       = cms.bool(True),
    doLaserCalFilt = cms.bool(False),

    jet_cuts                 = cms.bool(False),
    jet_minpt                = cms.double(30.0),
    jet_maxeta               = cms.double(2.4),
    min_jet                  = cms.int32(0),
    max_jet                  = cms.int32(4000),
    leading_jet_pt           = cms.double(30.0),

    muon_cuts                = cms.bool(True),
    tight_muon_minpt         = cms.double(10.0), # 26 GeV
    tight_muon_maxeta        = cms.double(2.1),
    loose_muon_minpt         = cms.double(10.0),
    loose_muon_maxeta        = cms.double(2.4),
    min_tight_muon           = cms.int32(0),

    electron_cuts                = cms.bool(True),
    tight_electron_minpt         = cms.double(0.0), # 30 GeV
    tight_electron_maxeta        = cms.double(999.),
    loose_electron_minpt         = cms.double(9999.0),
    loose_electron_maxeta        = cms.double(2.5),
    min_tight_electron           = cms.int32(0),
    
    min_tight_lepton         = cms.int32(0),
    max_tight_lepton         = cms.int32(999),
    trigger_consistent       = cms.bool(False),
    second_lepton_veto       = cms.bool(False),
    
    met_cuts                 = cms.bool(True),
    min_met                  = cms.double(0.0),
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
    electron_collection      = cms.InputTag('selectedPatElectronsPFlowLoose'),
    met_collection           = cms.InputTag('patMETsPFlow'),
    type1corrmet_collection  = cms.InputTag('pfType1CorrectedMet'),

    )


#######################################################
#
# Input files
#
process.inputs = cms.PSet (
   nEvents    = cms.int32(100000),
    skipEvents = cms.int32(0),
    lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange()),
    fileNames  = cms.vstring(
                           'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingletopWprime_M-1000_right_TuneZ2star_8TeV-comphep/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/SingletopWprime_M-1000_right_TuneZ2star_8TeV-comphep/USER/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/0000/FCB6CDEA-1F5E-E211-AEFA-003048FFD7D4.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingletopWprime_M-1000_right_TuneZ2star_8TeV-comphep/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/SingletopWprime_M-1000_right_TuneZ2star_8TeV-comphep/USER/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/0000/F86C3129-215E-E211-9E00-003048FFCBA4.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingletopWprime_M-1000_right_TuneZ2star_8TeV-comphep/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/SingletopWprime_M-1000_right_TuneZ2star_8TeV-comphep/USER/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/0000/F0F2D01A-215E-E211-89D6-003048FFD752.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingletopWprime_M-1000_right_TuneZ2star_8TeV-comphep/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/SingletopWprime_M-1000_right_TuneZ2star_8TeV-comphep/USER/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/0000/F0C16714-215E-E211-A020-003048FFD7A2.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingletopWprime_M-1000_right_TuneZ2star_8TeV-comphep/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/SingletopWprime_M-1000_right_TuneZ2star_8TeV-comphep/USER/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/0000/EC0F0929-215E-E211-8815-003048FFCBA4.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingletopWprime_M-1000_right_TuneZ2star_8TeV-comphep/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/SingletopWprime_M-1000_right_TuneZ2star_8TeV-comphep/USER/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/0000/D8532E29-215E-E211-BEBB-003048FFCBA4.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingletopWprime_M-1000_right_TuneZ2star_8TeV-comphep/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/SingletopWprime_M-1000_right_TuneZ2star_8TeV-comphep/USER/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/0000/D28C0830-215E-E211-94E0-003048FFCB96.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingletopWprime_M-1000_right_TuneZ2star_8TeV-comphep/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/SingletopWprime_M-1000_right_TuneZ2star_8TeV-comphep/USER/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/0000/CA2D7110-215E-E211-870D-003048FFD796.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingletopWprime_M-1000_right_TuneZ2star_8TeV-comphep/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/SingletopWprime_M-1000_right_TuneZ2star_8TeV-comphep/USER/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/0000/C62722F7-1F5E-E211-9227-003048FF9AA6.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingletopWprime_M-1000_right_TuneZ2star_8TeV-comphep/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/SingletopWprime_M-1000_right_TuneZ2star_8TeV-comphep/USER/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/0000/C26B5EE9-1F5E-E211-9F01-003048FFCB84.root',

     )
)


# JSON
JsonFile = ''
myList   = LumiList.LumiList(filename=JsonFile).getCMSSWString().split(',')
if (not process.ljmet.isMc == cms.bool(True)):
    process.inputs.lumisToProcess.extend(myList)
        
#######################################################
#
# Output
#
process.outputs = cms.PSet (
    outputName = cms.string('Wprime1000Right_noCuts_noTrig_noId_PFlowLoose_el'),
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
process.pfMuonSelector.cutsToIgnore = cms.vstring("maxPfRelIso")

# Loose muon                                                       
process.looseMuonSelector = process.pfMuonSelector.clone()
process.looseMuonSelector.version = cms.string('TOPPAG12_LJETS_VETO')

# electron
process.load('LJMet.Com.cutbasedIDSelector_cfi')
process.cutbasedIDSelector.version = cms.string('TIGHT')
process.cutbasedIDSelector.cutsToIgnore = cms.vstring("deta_EB", "dphi_EB", "sihih_EB", "hoe_EB", "d0_EB", "dZ_EB", "ooemoop_EB", "reliso_EB", "deta_EE", "dphi_EE", "sihih_EE", "hoe_EE", "d0_EB", "dZ_EB", "ooemoop_EE", "reliso_EE", "mHits", "vtxFitConv")

# loose electron
process.looseElectronSelector = process.cutbasedIDSelector.clone()
process.looseElectronSelector.version = cms.string('TIGHT')
process.looseElectronSelector.cutsToIgnore = cms.vstring()

# jets
process.load('PhysicsTools.SelectorUtils.pfJetIDSelector_cfi') 
process.pfJetIDSelector.version = cms.string('FIRSTDATA')
process.pfJetIDSelector.quality = cms.string('LOOSE')
