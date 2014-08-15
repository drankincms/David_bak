import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.Types as CfgTypes

inputs = cms.PSet (
    nEvents    = cms.int32(100000000),
    skipEvents = cms.int32(0),
    lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange()),
    fileNames  = cms.vstring()
    )

inputs.fileNames.extend([
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_100_1_JvB.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_101_1_2ah.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_102_1_uUQ.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_103_1_bvI.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_104_1_vMs.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_105_1_xR2.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_106_1_qF5.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_107_1_FgZ.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_108_1_WVS.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_109_1_cca.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_10_1_Zp2.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_110_1_yh2.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_111_1_wH9.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_112_1_wIx.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/lpctlbsm/meloam/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2/c04f3b4fa74c8266c913b71e0c74901d/tlbsm_53x_v2_mc_113_1_rcD.root'
    
    ])
