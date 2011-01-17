import FWCore.ParameterSet.Config as cms

htFilter = cms.EDFilter("HtFilter", 
  src = cms.InputTag("selectedPatJetsPF"),
  minHT = cms.double(100.00),
  maxHT = cms.double(-1.0)
)

metSqrtHtFilter = cms.EDFilter("MetSqrtHtFilter", 
  src = cms.InputTag("selectedPatJetsPF"),
  metSrc = cms.InputTag("patMETsPF"),
  minCut = cms.double(13.5),
  maxCut = cms.double(-1.0)
)

rawHtFilter = cms.EDFilter("RawHtFilter", 
  src = cms.InputTag("selectedPatJetsPF"),
  minHT = cms.double(100.00),
  maxHT = cms.double(-1.0)
)
