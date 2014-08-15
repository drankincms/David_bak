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
    
    trigger_cut  = cms.bool(True),
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
    tight_electron_minpt         = cms.double(10.0), # 30 GeV
    tight_electron_maxeta        = cms.double(2.5),
    loose_electron_minpt         = cms.double(20.0),
    loose_electron_maxeta        = cms.double(2.5),
    min_tight_electron           = cms.int32(1),
    
    min_tight_lepton         = cms.int32(1),
    max_tight_lepton         = cms.int32(1),
    trigger_consistent       = cms.bool(True),
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
    electron_collection      = cms.InputTag('selectedPatElectronsPFlow'),
    met_collection           = cms.InputTag('patMETsPFlow'),
    type1corrmet_collection  = cms.InputTag('pfType1CorrectedMet'),

    )


#######################################################
#
# Input files
#
process.inputs = cms.PSet (
   nEvents    = cms.int32(50000),
    skipEvents = cms.int32(0),
    lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange()),
    fileNames  = cms.vstring(
                           'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingletopWprime_M-2500_right_TuneZ2star_8TeV-comphep/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/SingletopWprime_M-2500_right_TuneZ2star_8TeV-comphep/USER/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/0000/10B1A180-4C28-E211-A360-002618943977.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingletopWprime_M-2500_right_TuneZ2star_8TeV-comphep/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/SingletopWprime_M-2500_right_TuneZ2star_8TeV-comphep/USER/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/0000/2A7F2966-4C28-E211-A508-0030486790B8.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingletopWprime_M-2500_right_TuneZ2star_8TeV-comphep/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/SingletopWprime_M-2500_right_TuneZ2star_8TeV-comphep/USER/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/0000/3C580D8D-4C28-E211-982B-00304867BF18.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingletopWprime_M-2500_right_TuneZ2star_8TeV-comphep/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/SingletopWprime_M-2500_right_TuneZ2star_8TeV-comphep/USER/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/0000/263919AC-4C28-E211-8765-00304867BFBC.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingletopWprime_M-2500_right_TuneZ2star_8TeV-comphep/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/SingletopWprime_M-2500_right_TuneZ2star_8TeV-comphep/USER/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/0000/3AAEF186-4C28-E211-8737-0026189438F9.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingletopWprime_M-2500_right_TuneZ2star_8TeV-comphep/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/SingletopWprime_M-2500_right_TuneZ2star_8TeV-comphep/USER/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/0000/44584D2C-4C28-E211-A199-003048678B08.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingletopWprime_M-2500_right_TuneZ2star_8TeV-comphep/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/SingletopWprime_M-2500_right_TuneZ2star_8TeV-comphep/USER/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/0000/4A0F037D-4C28-E211-8B95-003048678C9A.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingletopWprime_M-2500_right_TuneZ2star_8TeV-comphep/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/SingletopWprime_M-2500_right_TuneZ2star_8TeV-comphep/USER/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/0000/4A498480-4C28-E211-997C-003048678B7C.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingletopWprime_M-2500_right_TuneZ2star_8TeV-comphep/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/SingletopWprime_M-2500_right_TuneZ2star_8TeV-comphep/USER/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/0000/4CD60566-4C28-E211-BFED-002618943924.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingletopWprime_M-2500_right_TuneZ2star_8TeV-comphep/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/SingletopWprime_M-2500_right_TuneZ2star_8TeV-comphep/USER/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/0000/608A2597-4C28-E211-AF57-003048678BE8.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingletopWprime_M-2500_right_TuneZ2star_8TeV-comphep/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/SingletopWprime_M-2500_right_TuneZ2star_8TeV-comphep/USER/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/0000/6413CA79-4C28-E211-967D-003048678BE6.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingletopWprime_M-2500_right_TuneZ2star_8TeV-comphep/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/SingletopWprime_M-2500_right_TuneZ2star_8TeV-comphep/USER/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/0000/68ADB69E-4C28-E211-844A-003048D15E02.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingletopWprime_M-2500_right_TuneZ2star_8TeV-comphep/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/SingletopWprime_M-2500_right_TuneZ2star_8TeV-comphep/USER/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/0000/8439E658-4C28-E211-A7C6-003048678B7C.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingletopWprime_M-2500_right_TuneZ2star_8TeV-comphep/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/SingletopWprime_M-2500_right_TuneZ2star_8TeV-comphep/USER/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/0000/6E440B35-4C28-E211-AED9-003048678F74.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingletopWprime_M-2500_right_TuneZ2star_8TeV-comphep/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/SingletopWprime_M-2500_right_TuneZ2star_8TeV-comphep/USER/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/0000/88D66137-4C28-E211-B434-0026189438CC.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingletopWprime_M-2500_right_TuneZ2star_8TeV-comphep/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/SingletopWprime_M-2500_right_TuneZ2star_8TeV-comphep/USER/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/0000/98D4BC46-4C28-E211-9E1A-00304867924E.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingletopWprime_M-2500_right_TuneZ2star_8TeV-comphep/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/SingletopWprime_M-2500_right_TuneZ2star_8TeV-comphep/USER/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/0000/F6115365-4C28-E211-996A-00304867920A.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingletopWprime_M-2500_right_TuneZ2star_8TeV-comphep/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/SingletopWprime_M-2500_right_TuneZ2star_8TeV-comphep/USER/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/0000/96120E87-4C28-E211-85E8-0026189438E8.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingletopWprime_M-2500_right_TuneZ2star_8TeV-comphep/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/SingletopWprime_M-2500_right_TuneZ2star_8TeV-comphep/USER/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/0000/98D3E82A-4C28-E211-B21F-002618943922.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingletopWprime_M-2500_right_TuneZ2star_8TeV-comphep/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/SingletopWprime_M-2500_right_TuneZ2star_8TeV-comphep/USER/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/0000/BC918B78-4C28-E211-9972-00304867C026.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingletopWprime_M-2500_right_TuneZ2star_8TeV-comphep/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/SingletopWprime_M-2500_right_TuneZ2star_8TeV-comphep/USER/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/0000/DCE0B43A-4C28-E211-ABAC-003048678B7C.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingletopWprime_M-2500_right_TuneZ2star_8TeV-comphep/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/SingletopWprime_M-2500_right_TuneZ2star_8TeV-comphep/USER/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/0000/F8C70E24-4C28-E211-BEE9-002618FDA263.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingletopWprime_M-2500_right_TuneZ2star_8TeV-comphep/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/SingletopWprime_M-2500_right_TuneZ2star_8TeV-comphep/USER/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/0000/E09DC266-4C28-E211-99B4-0030486792AC.root',
                 'dcap:///pnfs/cms/WAX/11/store/results/B2G/SingletopWprime_M-2500_right_TuneZ2star_8TeV-comphep/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/SingletopWprime_M-2500_right_TuneZ2star_8TeV-comphep/USER/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2-c04f3b4fa74c8266c913b71e0c74901d/0000/FEBCA780-4C28-E211-A0C6-003048678ADA.root',

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
    outputName = cms.string('Wprime2500Right_noCuts_el'),
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
process.cutbasedIDSelector.cutsToIgnore     = cms.vstring("reliso_EB","reliso_EE")

# loose electron
process.looseElectronSelector = process.cutbasedIDSelector.clone()
process.looseElectronSelector.version = cms.string('VETO')
process.looseElectronSelector.cutsToIgnore     = cms.vstring()

# jets
process.load('PhysicsTools.SelectorUtils.pfJetIDSelector_cfi') 
process.pfJetIDSelector.version = cms.string('FIRSTDATA')
process.pfJetIDSelector.quality = cms.string('LOOSE')
