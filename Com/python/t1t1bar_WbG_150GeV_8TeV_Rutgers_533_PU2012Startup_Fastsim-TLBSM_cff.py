import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.Types as CfgTypes

inputs = cms.PSet (
    nEvents    = cms.int32(100000000),
    skipEvents = cms.int32(0),
    lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange()),
    fileNames  = cms.vstring()
    )

inputs.fileNames.extend([
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_150GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_10_1_QoR.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_150GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_11_1_ge6.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_150GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_12_1_6vy.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_150GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_13_1_tRl.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_150GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_14_1_nCL.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_150GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_15_1_ySk.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_150GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_16_1_KFf.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_150GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_17_1_Q0P.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_150GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_18_1_Bnl.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_150GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_19_1_THc.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_150GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_1_1_D2s.root ',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_150GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_20_1_OO3.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_150GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_2_1_7OE.root ',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_150GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_3_1_i8v.root ',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_150GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_4_1_EGS.root ',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_150GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_5_1_ayu.root ',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_150GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_6_1_kb6.root ',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_150GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_7_1_PAS.root ',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_150GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_8_1_XO2.root ',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_150GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_9_1_JI5.root '
    ])
