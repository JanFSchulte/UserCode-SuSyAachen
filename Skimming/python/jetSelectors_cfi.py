import FWCore.ParameterSet.Config as cms

patJetCountFilter = cms.EDFilter("PATJetCountFilter",
  src = cms.InputTag("basicJets"),
  minNumber = cms.uint32(1),
  cut = cms.string(""),
)
