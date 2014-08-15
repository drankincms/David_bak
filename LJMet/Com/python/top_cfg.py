import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes

process = cms.Process("LJMetCom")



############################################################
#
# FWLite application options
#
process.ljmet = cms.PSet(
    isMc      = cms.bool(True),
    runs = cms.vint32([])
)



############################################################
#
# Event selector options
#
process.event_selector = cms.PSet(
    # input parameter sets
    selection      = cms.string('TopSelector'),
    version        = cms.string('SelV3'),
    pvSrc          = cms.InputTag('goodOfflinePrimaryVertices'),
    muonSrc        = cms.InputTag('selectedPatMuonsPFlow'),
    electronSrc    = cms.InputTag('selectedPatElectronsPFlow'),

    muonStringSource = cms.string('selectedPatMuonsPFlow'),
    electronStringSource = cms.string('selectedPatElectronsPFlow'),


    jetSrc         = cms.InputTag('goodPatJetsPFlow'),
    jetSrcPF       = cms.InputTag('goodPatJetsPFlow'),

    jetStringSource = cms.string('goodPatJetsPFlow'),


    metSrc         = cms.InputTag('patMETPFlow'),
    metTCSrc       = cms.InputTag('patMETsPFlow'),
    metPFSrc       = cms.InputTag('patMETsPFlow'),
    metPfTypeISrc  = cms.InputTag('patMETsPFlow'),
    metMuonSrc     = cms.InputTag('patMETsPFlow'),
    metTypeIISrc   = cms.InputTag('patMETsPFlow'),

    metStringSource = cms.string('patMETsPFlow'),
    
    #rhoCorrection = cms.InputTag('kt6PFJetsPFlow:rho'),
    rhoCorrection = cms.InputTag('kt6PFJets:rho'),

    #triggerSrc     = cms.InputTag('patTrigger'),
    triggerSrc     = cms.InputTag('HLT_IsoMu24_v16'),
    BeamspotSrc    = cms.InputTag('offlineBeamSpot'),
    jetBTagSrc     = cms.InputTag('jetBProbabilityBJetTags'),
    jetFlavourSrc  = cms.InputTag('AK5byValAlgo'),
    tagInfo        = cms.InputTag('impactParameterTagInfos'),
    jetPModuleName = cms.InputTag('jetProbabilityBJetTags'),
    GenParticleSrc = cms.InputTag('prunedGenParticles'),


    useJEC =  cms.bool(False),
    JECup  =  cms.bool(True),

    usePU = cms.bool(True),
    useBTagSF = cms.bool(False),
    useLeptonIso = cms.bool(True),
    useTriggEff = cms.bool(True),
    useWxsec = cms.bool(False),
    
    JERup =   cms.bool(False),
    JERdown =   cms.bool(False),

    
    useBTagUncert =   cms.bool(False),
    BTagUncertUp =   cms.bool(False),
    BTagUncertDown =   cms.bool(False),


    removeOverlapEvents = cms.bool(False)

    )




#######################################################
#
# below are object selectors' options
#


# legacy
# we don't use this selector for
# this top analysis but keep it
# for backwards compatibility
#process.load('LJMet.Com.pvSelector_cfi')
process.load('PhysicsTools.SelectorUtils.pvSelector_cfi')

# this is the actual PV selector that we use
process.pvObjectSelector = cms.PSet(
    version = cms.string('DATA2010'),
    quality = cms.string('Top_SelV3')
    )

process.muonSelector = cms.PSet(
    version = cms.string('DATA2010'),
    quality = cms.string('TOP_2012')
    )


process.looseMuonSelector = process.muonSelector.clone(
    quality = cms.string('Top_Loose_SelV3')
    )

process.electronSelector = cms.PSet(
    version = cms.string('DATA2010'),
    quality = cms.string('Top_SelV3')
    )    

process.looseElectronSelector = process.electronSelector.clone(
    version = cms.string('DATA2010'),
    quality = cms.string('Top_SelV3')
    )


process.jetIDSelector = cms.PSet(
    version = cms.string('DATA2010'),
    quality = cms.string('Top_SelV3')
    )



process.metSelector = cms.PSet(
    version = cms.string('DATA2010'),
    quality = cms.string('Off')
    )



###############################################################
#
# Source
#

process.inputs = cms.PSet (
    nEvents    = cms.int32(5000),
    skipEvents = cms.int32(0),
    fileNames  = cms.vstring(
        
        #'dcap:////pnfs/cms/WAX/11/store//user//lpctlbsm/guragain/TTJets_TuneZ2star_8TeV-madgraph-tauola/8TeV_Summer12-PU_S7_START52_V5-v1/e26ff2b873f6dfa53bf7c2f73f519ef0/ttbsm_52x_mc_17_1_Jcw.root'
        'file:/uscms_data/d2/mike1886/LJMet/CMSSW_5_2_5/src/TopQuarkAnalysis/TopPairBSM/test/ttbsm_52x_mc.root'

        ),
    lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
)


# get JSON file correctly parced
#JSONfile = '../bin/JSON_FILES/Cert_190456-194076_8TeV_PromptReco_Collisions12_JSON.txt'
JSONfile = 'data/json/Cert_160404-180252_7TeV_ReRecoNov08_Collisions11_JSON_MuonPhys_v2.txt'
myList = LumiList.LumiList (filename = JSONfile).getCMSSWString().split(',')

if not process.ljmet.isMc:
    process.inputs.lumisToProcess.extend(myList)

###############################################################
#
# Output
#

process.outputs = cms.PSet (
    outputName = cms.string('topPlotsTTbar.root'),
    treeName   = cms.string('treetop_mu'),
    treeName_el   = cms.string('treetop_el')
)
