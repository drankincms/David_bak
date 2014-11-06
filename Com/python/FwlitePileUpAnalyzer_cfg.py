import FWCore.ParameterSet.Config as cms
 
process = cms.Process("FWLitePlots")

process.PileUpStudies = cms.PSet(
    # input parameter sets
    pvSrc = cms.InputTag('offlinePrimaryVertices'),
    muonSrc = cms.InputTag('selectedPatMuons'),
    metSrc = cms.InputTag('patMETs'),   
    jetSrc = cms.InputTag('selectedPatJets'),
    jetBTagSrc =  cms.InputTag('jetBProbabilityBJetTags'),
    jetFlavourSrc =  cms.InputTag('AK5byValAlgo')
    )



process.load('LJMet.Com.muonSelector_cfi')
process.load('LJMet.Com.jetIDSelector_cfi')
process.load('LJMet.Com.metSelector_cfi')
process.load('LJMet.Com.pvSelector_cfi')


process.plotParameters = cms.PSet (
    useMC = cms.bool(False),
    runs = cms.vint32([])
)


import FWCore.Utilities.FileUtils as FileUtils


process.inputs = cms.PSet (
    fileNames = cms.vstring(

    #'dcap:///pnfs/fnal.gov/usr/cms/WAX/resilient/mike1886/PileUpStudies/W_Jets_NoPileUp/PatTuple_1_1_Pbt.root',
    #'dcap:///pnfs/fnal.gov/usr/cms/WAX/resilient/mike1886/PileUpStudies/W_Jets_NoPileUp/PatTuple_10_1_cTK.root',
    #'dcap:///pnfs/fnal.gov/usr/cms/WAX/resilient/mike1886/PileUpStudies/W_Jets_NoPileUp/PatTuple_100_1_Qu6.root',

    'dcap:///pnfs/fnal.gov/usr/cms/WAX/resilient/mike1886/PileUpStudies/W_Jets_WithPileUp/PatTuple_1_1_eOG.root',
    'dcap:///pnfs/fnal.gov/usr/cms/WAX/resilient/mike1886/PileUpStudies/W_Jets_WithPileUp/PatTuple_10_1_1gM.root',
    'dcap:///pnfs/fnal.gov/usr/cms/WAX/resilient/mike1886/PileUpStudies/W_Jets_WithPileUp/PatTuple_100_1_wyM.root',


    )
)
 
process.outputs = cms.PSet (
    
    outputName = cms.string('./WWithPileup.root')
  
)

process.maxEvents = cms.PSet (
    number = cms.int32(1000)
)

