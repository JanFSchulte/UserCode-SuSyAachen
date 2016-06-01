import FWCore.ParameterSet.Config as cms

mhtFilter = cms.EDFilter("MhtFilter", 
  src = cms.InputTag("cleanPatJetsAK4Calo"),
  minMHT = cms.double(0.0)
)
