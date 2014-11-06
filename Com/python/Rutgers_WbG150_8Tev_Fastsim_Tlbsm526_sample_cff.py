import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.Types as CfgTypes

inputs = cms.PSet (
    nEvents    = cms.int32(100000000),
    skipEvents = cms.int32(0),
    lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange()),
    fileNames  = cms.vstring()
    )

inputs.fileNames.extend([
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_150GeV_8TeV_Rutgers_nocuts_v1_Fastsim-AODSIM/TLBSM_v1-PAT/ttbsm_52x_mc_10_1_JjT.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_150GeV_8TeV_Rutgers_nocuts_v1_Fastsim-AODSIM/TLBSM_v1-PAT/ttbsm_52x_mc_11_1_Na7.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_150GeV_8TeV_Rutgers_nocuts_v1_Fastsim-AODSIM/TLBSM_v1-PAT/ttbsm_52x_mc_12_1_h7k.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_150GeV_8TeV_Rutgers_nocuts_v1_Fastsim-AODSIM/TLBSM_v1-PAT/ttbsm_52x_mc_13_1_D24.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_150GeV_8TeV_Rutgers_nocuts_v1_Fastsim-AODSIM/TLBSM_v1-PAT/ttbsm_52x_mc_14_1_7W6.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_150GeV_8TeV_Rutgers_nocuts_v1_Fastsim-AODSIM/TLBSM_v1-PAT/ttbsm_52x_mc_15_1_wcd.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_150GeV_8TeV_Rutgers_nocuts_v1_Fastsim-AODSIM/TLBSM_v1-PAT/ttbsm_52x_mc_16_1_L0Y.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_150GeV_8TeV_Rutgers_nocuts_v1_Fastsim-AODSIM/TLBSM_v1-PAT/ttbsm_52x_mc_17_1_t65.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_150GeV_8TeV_Rutgers_nocuts_v1_Fastsim-AODSIM/TLBSM_v1-PAT/ttbsm_52x_mc_18_1_tKs.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_150GeV_8TeV_Rutgers_nocuts_v1_Fastsim-AODSIM/TLBSM_v1-PAT/ttbsm_52x_mc_19_1_WOI.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_150GeV_8TeV_Rutgers_nocuts_v1_Fastsim-AODSIM/TLBSM_v1-PAT/ttbsm_52x_mc_1_1_tnQ.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_150GeV_8TeV_Rutgers_nocuts_v1_Fastsim-AODSIM/TLBSM_v1-PAT/ttbsm_52x_mc_20_1_IDN.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_150GeV_8TeV_Rutgers_nocuts_v1_Fastsim-AODSIM/TLBSM_v1-PAT/ttbsm_52x_mc_2_1_Hjz.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_150GeV_8TeV_Rutgers_nocuts_v1_Fastsim-AODSIM/TLBSM_v1-PAT/ttbsm_52x_mc_3_1_4kY.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_150GeV_8TeV_Rutgers_nocuts_v1_Fastsim-AODSIM/TLBSM_v1-PAT/ttbsm_52x_mc_4_1_uu8.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_150GeV_8TeV_Rutgers_nocuts_v1_Fastsim-AODSIM/TLBSM_v1-PAT/ttbsm_52x_mc_5_1_uLr.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_150GeV_8TeV_Rutgers_nocuts_v1_Fastsim-AODSIM/TLBSM_v1-PAT/ttbsm_52x_mc_6_1_fsv.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_150GeV_8TeV_Rutgers_nocuts_v1_Fastsim-AODSIM/TLBSM_v1-PAT/ttbsm_52x_mc_7_1_5GV.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_150GeV_8TeV_Rutgers_nocuts_v1_Fastsim-AODSIM/TLBSM_v1-PAT/ttbsm_52x_mc_8_1_8dt.root',
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kukarzev/t1t1bar_WbG_150GeV_8TeV_Rutgers_nocuts_v1_Fastsim-AODSIM/TLBSM_v1-PAT/ttbsm_52x_mc_9_1_HaZ.root'
    ])
