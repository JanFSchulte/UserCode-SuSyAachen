import FWCore.ParameterSet.Config as cms

defaultSelector = cms.EDFilter("PATIsoTrackSelector", 
    filter = cms.bool(True),
    src = cms.InputTag("isolatedTracks"),
    name = cms.InputTag("basicIsolatedTracks"),
    muonSource = cms.InputTag("slimmedMuons"),
    electronSource = cms.InputTag("slimmedElectrons"),
    pfCandSource = cms.InputTag("packedPFCandidates")
    )
