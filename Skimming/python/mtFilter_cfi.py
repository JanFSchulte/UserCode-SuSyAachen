import FWCore.ParameterSet.Config as cms

mtFilter = cms.EDFilter("MtFilter", 
  src = cms.InputTag("selectedPatMuonsPF"),
  srcMET = cms.InputTag("patMETsPF"),
  minMT = cms.double(30.00),
  maxMT = cms.double(-1)
)
