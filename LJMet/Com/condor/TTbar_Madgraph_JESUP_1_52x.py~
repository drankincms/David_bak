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
process.WprimeCalc.isTB = cms.bool(False)

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
    min_tight_muon           = cms.int32(0),

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
    JECup                    = cms.bool(True),
    JECdown                  = cms.bool(False),
    JERup                    = cms.bool(False),
    JERdown                  = cms.bool(False),
    JEC_txtfile = cms.string('../cond/Fall12_V6_DATA_UncertaintySources_AK5PFchs.txt'),
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
   nEvents    = cms.int32(200000),
    skipEvents = cms.int32(0),
    lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange()),
    fileNames  = cms.vstring(
                           'dcap:///pnfs/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_1_1_v0T.root',
                 'dcap:///pnfs/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_2_1_y3T.root',
                 'dcap:///pnfs/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_5_1_Qpb.root',
                 'dcap:///pnfs/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_101_1_2ah.root',
                 'dcap:///pnfs/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_100_1_JvB.root',
                 'dcap:///pnfs/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_103_1_bvI.root',
                 'dcap:///pnfs/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_102_1_uUQ.root',
                 'dcap:///pnfs/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_105_1_xR2.root',
                 'dcap:///pnfs/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_104_1_vMs.root',
                 'dcap:///pnfs/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_107_1_FgZ.root',
                 'dcap:///pnfs/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_106_1_qF5.root',
                 'dcap:///pnfs/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_109_1_cca.root',
                 'dcap:///pnfs/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_108_1_WVS.root',
                 'dcap:///pnfs/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_110_1_yh2.root',
                 'dcap:///pnfs/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_10_1_Zp2.root',
                 'dcap:///pnfs/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_112_1_wIx.root',
                 'dcap:///pnfs/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_111_1_wH9.root',
                 'dcap:///pnfs/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_114_1_T1C.root',
                 'dcap:///pnfs/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_113_1_rcD.root',
                 'dcap:///pnfs/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_116_1_X2Z.root',
                 'dcap:///pnfs/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_115_1_rGV.root',
                 'dcap:///pnfs/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_118_1_B9m.root',
                 'dcap:///pnfs/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_117_1_S1c.root',
                 'dcap:///pnfs/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_11_1_DSp.root',
                 'dcap:///pnfs/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_119_1_lDH.root',
                 'dcap:///pnfs/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_121_1_rbq.root',
                 'dcap:///pnfs/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_120_1_9MY.root',
                 'dcap:///pnfs/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_123_1_W7Z.root',
                 'dcap:///pnfs/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_122_1_8Fn.root',
                 'dcap:///pnfs/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_125_1_73x.root',
                 'dcap:///pnfs/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_124_1_WSK.root',
                 'dcap:///pnfs/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_127_1_prd.root',
                 'dcap:///pnfs/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_126_1_TC8.root',
                 'dcap:///pnfs/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_129_1_MOK.root',
                 'dcap:///pnfs/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_128_1_RaT.root',
                 'dcap:///pnfs/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_130_1_gMZ.root',
                 'dcap:///pnfs/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_12_2_cxO.root',
                 'dcap:///pnfs/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_132_1_Utx.root',
                 'dcap:///pnfs/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_131_1_DA9.root',
                 'dcap:///pnfs/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_134_1_B9K.root',
                 'dcap:///pnfs/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_133_1_htO.root',
                 'dcap:///pnfs/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_136_1_E3J.root',
                 'dcap:///pnfs/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_135_1_jZC.root',
                 'dcap:///pnfs/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_138_1_yKL.root',
                 'dcap:///pnfs/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_137_1_VgP.root',
                 'dcap:///pnfs/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_13_1_Vxu.root',
                 'dcap:///pnfs/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_139_1_szn.root',
                 'dcap:///pnfs/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_141_1_RhD.root',
                 'dcap:///pnfs/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_140_1_tVL.root',
                 'dcap:///pnfs/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_143_1_kIh.root',

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
    outputName = cms.string('TTbar_Madgraph_JESUP_1'),
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
