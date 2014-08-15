import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.Types as CfgTypes

inputs = cms.PSet (
    nEvents    = cms.int32(100000000),
    skipEvents = cms.int32(0),
    lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange()),
    fileNames  = cms.vstring()
    )

inputs.fileNames.extend([
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_120GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_10_1_Mgh.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_120GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_11_1_MDe.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_120GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_12_1_F2d.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_120GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_13_1_Ajb.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_120GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_14_1_fHF.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_120GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_15_1_Vwh.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_120GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_16_1_0CM.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_120GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_17_1_Q4a.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_120GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_18_1_UTA.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_120GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_19_1_uMa.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_120GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_1_1_THf.root ',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_120GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_20_1_SGW.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_120GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_2_1_lA3.root ',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_120GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_3_1_p9R.root ',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_120GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_4_1_uqB.root ',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_120GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_5_1_Eip.root ',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_120GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_6_1_Sma.root ',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_120GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_7_1_x3k.root ',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_120GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_8_1_YrM.root ',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_120GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_9_1_OVv.root '
    ])
