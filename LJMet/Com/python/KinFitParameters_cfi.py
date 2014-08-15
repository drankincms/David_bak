import FWCore.ParameterSet.Config as cms

defaultKinFitParameters = cms.PSet(
    kinFitterLabel = cms.string("kinFitSemiLepEventCSVm"),
    #    jets = cms.InputTag("selectedPatJets"),
    #    leps = cms.InputTag("selectedPatElectrons"),
    #    mets = cms.InputTag("patMETs"),

    # ------------------------------------------------
    # maximum number of jets to be considered in the
    # jet combinatorics (has to be >= 4, can be set to
    # -1 if you want to take all)
    # ------------------------------------------------
    maxNJets = cms.int32(5),
    
    #-------------------------------------------------
    # maximum number of jet combinations finally
    # written into the event, starting from the "best"
    # (has to be >= 1, can be set to -1 if you want to
    # take all)
    #-------------------------------------------------
    maxNComb = cms.int32(-1),
    
    # ------------------------------------------------
    # option to take only a given jet combination
    # instead of going through the full combinatorics
    # ------------------------------------------------
    match = cms.InputTag(""),
    useOnlyMatch = cms.bool(False),
    
    # ------------------------------------------------
    # option to use b-tagging
    # ------------------------------------------------
    #bTagAlgo          = cms.string("trackCountingHighEffBJetTags"),
    #minBDiscBJets     = cms.double(3.3),
    #maxBDiscLightJets = cms.double(3.3),
    bTagAlgo          = cms.string("combinedSecondaryVertexBJetTags"),
    minBDiscBJets     = cms.double(0.679),
    maxBDiscLightJets = cms.double(0.679),
    useBTagging       = cms.bool(True),
    
    # ------------------------------------------------
    # settings for the KinFitter
    # ------------------------------------------------
    maxNrIter = cms.uint32(500),
    maxDeltaS = cms.double(5e-05),
    maxF      = cms.double(0.0001),
    # ------------------------------------------------
    # select parametrisation
    # 0: EMom, 1: EtEtaPhi, 2: EtThetaPhi
    # ------------------------------------------------
    jetParametrisation = cms.uint32(1),
    lepParametrisation = cms.uint32(1),
    metParametrisation = cms.uint32(1),

    # ------------------------------------------------
    # set constraints
    # 1: Whad-mass, 2: Wlep-mass, 3: thad-mass,
    # 4: tlep-mass, 5: nu-mass, 6: equal t-masses
    # ------------------------------------------------
    constraints = cms.vuint32(1, 2, 6),
    
    # ------------------------------------------------
    # set mass values used in the constraints
    # ------------------------------------------------
    mW   = cms.double(80.4),
    mTop = cms.double(173.)
    )
