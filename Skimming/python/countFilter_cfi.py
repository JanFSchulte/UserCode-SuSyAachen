import FWCore.ParameterSet.Config as cms

candCountFilter = cms.EDFilter("CandCountFilter", 
  src = cms.InputTag("selectedPatJetsPF"),
  minNumber = cms.uint32(0),
  maxNumber = cms.uint32(9999999)
)

vertexCountFilter = cms.EDFilter("nVertexCountFilter", 
  src = cms.InputTag("selectedPatJetsPF"),
  minNumber = cms.uint32(0),
  maxNumber = cms.uint32(9999999)
)

candViewCountFilter = cms.EDFilter("CandViewCountFilter",
  src = cms.InputTag("selectedPatJetsPF"),
  minNumber = cms.uint32(0),
  maxNumber = cms.uint32(9999999)
)
