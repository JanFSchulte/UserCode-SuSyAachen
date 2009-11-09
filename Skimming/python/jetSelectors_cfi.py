import FWCore.ParameterSet.Config as cms

patJetCountFilter = cms.EDFilter("PATJetCountFilter", filter = cms.bool(False),
  src = cms.InputTag("basicJets"),
  minNumber = cms.uint32(1),
  cut = cms.string(""),
)
