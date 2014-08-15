import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.Types as CfgTypes

inputs = cms.PSet (
    nEvents    = cms.int32(100000000),
    skipEvents = cms.int32(0),
    lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange()),
    fileNames  = cms.vstring()
    )

inputs.fileNames.extend([
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_180GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_10_1_K3z.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_180GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_11_1_zFL.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_180GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_12_1_mJW.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_180GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_13_1_YcS.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_180GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_14_1_mao.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_180GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_15_1_A03.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_180GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_16_2_ErG.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_180GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_17_2_SD0.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_180GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_18_1_b0g.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_180GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_19_1_kPg.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_180GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_1_1_sgl.root ',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_180GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_20_1_W7x.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_180GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_2_1_Xle.root ',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_180GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_3_1_18J.root ',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_180GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_4_1_JUT.root ',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_180GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_5_1_dQg.root ',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_180GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_6_1_FF5.root ',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_180GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_7_1_cn9.root ',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_180GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_8_1_OJ3.root ',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_180GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-AODSIM/TLBSM_53X_v2/tlbsm_53x_v2_mc_9_1_ODt.root '
    ])
