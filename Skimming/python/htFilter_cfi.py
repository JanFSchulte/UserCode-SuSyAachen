import FWCore.ParameterSet.Config as cms

htFilter = cms.EDFilter("HtFilter", 
  src = cms.InputTag("selectedPatJetsPF"),
  minHT = cms.double(100.00),
  maxHT = cms.double(-1.0)
)
