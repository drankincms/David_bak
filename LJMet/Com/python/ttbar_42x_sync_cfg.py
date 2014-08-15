import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes

process = cms.Process("ttbar")

_useMC = True



############################################################
#
# Event selector options
#
process.event_selector = cms.PSet(

    selection = cms.string('ttbar'),

    # cuts
    trigger_cut  = cms.bool(True),
    dump_trigger = cms.bool(True),
    trigger_path = cms.string('HLT_IsoMu17_v5'),
    
    pv_cut         = cms.bool(True),

    hbhe_cut       = cms.bool(False),

    jet_cuts                 = cms.bool(True),
    jet_minpt                = cms.double(30.0),
    jet_maxeta               = cms.double(2.4),
    min_jet                  = cms.int32(4),
    max_jet                  = cms.int32(4000),

    muon_cuts                = cms.bool(True),
    min_tight_muon           = cms.int32(1),
    tight_muon_minpt         = cms.double(20.0),
    tight_muon_maxeta        = cms.double(2.1),
    tight_muon_mindeltaR_jet = cms.double(0.3),
    max_tight_muon           = cms.int32(1),

    loose_muon_minpt         = cms.double(10.0),
    loose_muon_maxeta        = cms.double(2.5),
    loose_muon_veto          = cms.bool(True),

    electron_minEt           = cms.double(15.0),
    electron_maxeta          = cms.double(2.5),
    electron_veto            = cms.bool(True),

    met_cuts                 = cms.bool(False),
    min_met                  = cms.double(0.0),

    btag_cuts                = cms.bool(True),
    btagger                  = cms.string('simpleSecondaryVertexHighEffBJetTags'),
    btag_min_discr           = cms.double(2.0),
    btag_1                   = cms.bool(True),
    btag_2                   = cms.bool(True),
    btag_3                   = cms.bool(False),

    # input collections
    trigger_collection       = cms.InputTag('TriggerResults::HLT'),
    pv_collection            = cms.InputTag('goodOfflinePrimaryVertices'),
    jet_collection           = cms.InputTag('goodPatJetsPFlow'),
    muon_collection          = cms.InputTag('selectedPatMuonsPFlow'),
    electron_collection      = cms.InputTag('selectedPatElectronsPFlow'),
    met_collection           = cms.InputTag('patMETsPFlow'),

    )



#######################################################
#
# Object selector options
#



# Primary vertex
process.load('PhysicsTools.SelectorUtils.pvSelector_cfi')
process.pvSelector.pvSrc   = cms.InputTag('goodOfflinePrimaryVertices')
process.pvSelector.minNdof = cms.double(4.0)
process.pvSelector.maxZ    = cms.double(24.0)
process.pvSelector.maxRho  = cms.double(2.0)



# tight muon
process.load('PhysicsTools.SelectorUtils.pfMuonSelector_cfi') 
process.pfMuonSelector.version          = cms.string('SPRING11')
process.pfMuonSelector.Chi2             = cms.double(10.0)
process.pfMuonSelector.NHits            = cms.int32(11)
process.pfMuonSelector.NValMuHits       = cms.int32(1)
process.pfMuonSelector.D0               = cms.double(0.02)
process.pfMuonSelector.PFIso            = cms.double(0.125)
process.pfMuonSelector.nPixelHits       = cms.int32(1)
process.pfMuonSelector.nMatchedStations = cms.int32(2)
process.pfMuonSelector.cutsToIgnore     = cms.vstring()



# loose muon
process.looseMuonSelector = process.pfMuonSelector.clone()
process.looseMuonSelector.PFIso        = cms.double(0.2)
process.looseMuonSelector.cutsToIgnore = cms.vstring('TrackerMuon',
                                                     'Chi2',
                                                     'NHits',
                                                     'NValMuHits',
                                                     'D0',
                                                     'nPixelHits',
                                                     'nMatchedStations')



# electron
process.load('PhysicsTools.SelectorUtils.pfElectronSelector_cfi')
process.pfElectronSelector.version = cms.string('SPRING11')
process.pfElectronSelector.PFIso = cms.double(0.2)
process.pfElectronSelector.cutsToIgnore = cms.vstring('electronID',
                                                      'MVA',
                                                      'MaxMissingHits',
                                                      'D0',
                                                      'electronIDused',
                                                      'ConversionRejection')



# jets
process.load('PhysicsTools.SelectorUtils.pfJetIDSelector_cfi') 
process.pfJetIDSelector.version = cms.string('FIRSTDATA')
process.pfJetIDSelector.quality = cms.string('LOOSE')







############################################################
#
# FWLite application options
#



# plotting
process.plotParameters = cms.PSet (
    useMC = cms.bool(_useMC),
    runs = cms.vint32([])
    )



# Source
process.inputs = cms.PSet (
    nEvents    = cms.int32(1000000),
    skipEvents = cms.int32(0),
    fileNames  = cms.vstring(
        # sync 2011 ttbar mu+jets
        'file:/eos/uscms/store/user/kukarzev/data/SevenTev/sync/top2011/ttbar_ttbsm_42x_mc_v3.root'
        ),
    lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
)       
        
        
        
# Output
process.outputs = cms.PSet (
    outputName = cms.string('topPlotsTTbar.root'),
    treeName   = cms.string('treetop_mu'),
    treeName_el   = cms.string('treetop_el')
)
