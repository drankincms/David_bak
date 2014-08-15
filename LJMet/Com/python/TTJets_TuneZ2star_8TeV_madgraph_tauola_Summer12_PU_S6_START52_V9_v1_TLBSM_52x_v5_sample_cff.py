import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.Types as CfgTypes

inputs = cms.PSet (
    nEvents    = cms.int32(100000000),
    skipEvents = cms.int32(0),
    lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange()),
    fileNames  = cms.vstring()
    )

inputs.fileNames.extend([
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/ra2tau/meloam/TTJets_TuneZ2star_8TeV-madgraph-tauola/Summer12-PU_S6_START52_V9-v1_TLBSM_52x_v5/2bcc6fdd1e664d93e9026c3764d0b403/ttbsm_52x_mc_9_1_rKE.root'
    ])
