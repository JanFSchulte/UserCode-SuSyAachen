import FWCore.ParameterSet.Config as cms

defaultSelector = cms.EDFilter("BadGlobalMuonTagger",
    src  = cms.InputTag("muons"),
    vtx   = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muonPtCut = cms.double(20),
    selectClones = cms.bool(True),
    taggingMode = cms.bool(False),
    verbose = cms.bool(False),
)


