import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")
process.load("PhysicsTools.HepMCCandAlgos.genParticles_cfi")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('SMSvalidation.root')
)



process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.cerr.threshold = 'INFO'
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.load("PhysicsTools.HepMCCandAlgos.genParticles_cfi")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    generator = cms.PSet(
        initialSeed = cms.untracked.uint32(123456789),
        engineName = cms.untracked.string('HepJamesRandom')
    )
)


process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    generator = cms.PSet(
        initialSeed = cms.untracked.uint32(123456789),
        engineName = cms.untracked.string('HepJamesRandom')
    )
)
#process.randomEngineStateProducer = cms.EDProducer("RandomEngineStateProducer")

process.source = cms.Source("LHESource",
    fileNames = cms.untracked.vstring('file:test.lhe')
)



process.generator = cms.EDFilter("Pythia6HadronizerFilter",
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    maxEventsToPrint = cms.untracked.int32(1),
    pythiaPylistVerbosity = cms.untracked.int32(1),
    comEnergy = cms.double(7000.0),
    PythiaParameters = cms.PSet(
        processParameters = cms.vstring('MSEL=0         ! User defined processes','MSTP(81)=0','MSTP(91)=0','MSTP(111)=0',
										'MSTP(61)=0','MSTP(71)=0',
                        'PMAS(5,1)=4.4   ! b quark mass',
                        'PMAS(6,1)=172.4 ! t quark mass',
                        'MSTJ(1)=       ! Fragmentation/hadronization on or off',
                        'MSTP(61)=      ! Parton showering on or off'),
        # This is a vector of ParameterSet names to be read, in this order
        parameterSets = cms.vstring(
            'processParameters')
    ),
)



process.demo = cms.EDAnalyzer('SimpleLHEAnalysis'
)


process.p = cms.Path(process.generator*process.genParticles*process.demo)
#process.p1 = cms.Path(process.randomEngineStateProducer)
#process.outpath = cms.EndPath(process.GEN)

#process.schedule = cms.Schedule(process.p, process.p1, process.outpath)
process.schedule = cms.Schedule(process.p)

