import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.Types as CfgTypes

inputs = cms.PSet (
    nEvents    = cms.int32(10000),
    skipEvents = cms.int32(0),
    lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange()),
    fileNames  = cms.vstring()
    )

inputs.fileNames.extend([

    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/results/B2G/SingleMu/StoreResults-SingleMu_Run2012A-13Jul2012-v1_TLBSM_53x_v2_jsonfix-e3fb55b810dc7a0811f4c66dfa2267c9/SingleMu/USER/StoreResults-SingleMu_Run2012A-13Jul2012-v1_TLBSM_53x_v2_jsonfix-e3fb55b810dc7a0811f4c66dfa2267c9/0000/00329188-173B-E211-9560-003048FFCB8C.root'


    #'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/lpctlbsm/knash/SingleMu/SingleMu_Run2012A-13Jul2012-v1_TLBSM_53x_v2/e3fb55b810dc7a0811f4c66dfa2267c9/tlbsm_53x_v2_data_9_1_3z8.root'
    ])
