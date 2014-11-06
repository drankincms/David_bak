import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.Types as CfgTypes

inputs = cms.PSet (
    nEvents    = cms.int32(100000000),
    skipEvents = cms.int32(0),
    lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange()),
    fileNames  = cms.vstring()
    )

inputs.fileNames.extend([
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_172GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_10_1_dXZ.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_172GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_11_1_7ZN.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_172GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_12_1_ODn.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_172GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_13_1_yOF.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_172GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_14_1_yNt.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_172GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_15_1_xrE.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_172GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_16_1_iWr.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_172GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_17_1_tFP.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_172GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_18_1_ggq.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_172GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_19_1_Uqn.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_172GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_1_1_vtz.root ',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_172GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_20_1_vvG.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_172GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_2_1_iaG.root ',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_172GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_3_1_zEG.root ',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_172GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_4_1_mGG.root ',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_172GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_5_1_nbt.root ',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_172GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_6_1_7Vw.root ',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_172GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_7_1_u4x.root ',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_172GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_8_1_1Zv.root ',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_172GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_9_1_ilf.root '
    ])
