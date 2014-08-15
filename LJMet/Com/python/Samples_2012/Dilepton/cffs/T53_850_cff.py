import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.Types as CfgTypes

inputs = cms.PSet (
    nEvents    = cms.int32(100000000),
    skipEvents = cms.int32(0),
    lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange()),
    fileNames  = cms.vstring()
    )
inputs.fileNames.extend([

    'dcache:/pnfs/cms/WAX/11//store/user/lpctlbsm/avetisya/T53/T53T53To2L2Nu_850_8TeV_Summer12/avetisya/T53T53To2L2Nu_M-850_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2_T53T53To2L2Nu_850_8TeV_Summer12/3abcc3b1cd74b7b45c7ed2df0ee1e03c//tlbsm_53x_v2_mc_1_1_hsb.root',
    'dcache:/pnfs/cms/WAX/11//store/user/lpctlbsm/avetisya/T53/T53T53To2L2Nu_850_8TeV_Summer12/avetisya/T53T53To2L2Nu_M-850_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2_T53T53To2L2Nu_850_8TeV_Summer12/3abcc3b1cd74b7b45c7ed2df0ee1e03c//tlbsm_53x_v2_mc_2_1_aD9.root',
    'dcache:/pnfs/cms/WAX/11//store/user/lpctlbsm/avetisya/T53/T53T53To2L2Nu_850_8TeV_Summer12/avetisya/T53T53To2L2Nu_M-850_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2_T53T53To2L2Nu_850_8TeV_Summer12/3abcc3b1cd74b7b45c7ed2df0ee1e03c//tlbsm_53x_v2_mc_3_1_POT.root',
    'dcache:/pnfs/cms/WAX/11//store/user/lpctlbsm/avetisya/T53/T53T53To2L2Nu_850_8TeV_Summer12/avetisya/T53T53To2L2Nu_M-850_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2_T53T53To2L2Nu_850_8TeV_Summer12/3abcc3b1cd74b7b45c7ed2df0ee1e03c//tlbsm_53x_v2_mc_4_1_pX1.root',
    'dcache:/pnfs/cms/WAX/11//store/user/lpctlbsm/avetisya/T53/T53T53To2L2Nu_850_8TeV_Summer12/avetisya/T53T53To2L2Nu_M-850_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2_T53T53To2L2Nu_850_8TeV_Summer12/3abcc3b1cd74b7b45c7ed2df0ee1e03c//tlbsm_53x_v2_mc_5_1_mX5.root',
    'dcache:/pnfs/cms/WAX/11//store/user/lpctlbsm/avetisya/T53/T53T53To2L2Nu_850_8TeV_Summer12/avetisya/T53T53To2L2Nu_M-850_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2_T53T53To2L2Nu_850_8TeV_Summer12/3abcc3b1cd74b7b45c7ed2df0ee1e03c//tlbsm_53x_v2_mc_6_1_k57.root',
    'dcache:/pnfs/cms/WAX/11//store/user/lpctlbsm/avetisya/T53/T53T53To2L2Nu_850_8TeV_Summer12/avetisya/T53T53To2L2Nu_M-850_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2_T53T53To2L2Nu_850_8TeV_Summer12/3abcc3b1cd74b7b45c7ed2df0ee1e03c//tlbsm_53x_v2_mc_7_1_EzO.root',
    'dcache:/pnfs/cms/WAX/11//store/user/lpctlbsm/avetisya/T53/T53T53To2L2Nu_850_8TeV_Summer12/avetisya/T53T53To2L2Nu_M-850_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2_T53T53To2L2Nu_850_8TeV_Summer12/3abcc3b1cd74b7b45c7ed2df0ee1e03c//tlbsm_53x_v2_mc_8_1_mUr.root',
    'dcache:/pnfs/cms/WAX/11//store/user/lpctlbsm/avetisya/T53/T53T53To2L2Nu_850_8TeV_Summer12/avetisya/T53T53To2L2Nu_M-850_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v2_T53T53To2L2Nu_850_8TeV_Summer12/3abcc3b1cd74b7b45c7ed2df0ee1e03c//tlbsm_53x_v2_mc_9_1_pfF.root',    ])
