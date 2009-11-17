import FWCore.ParameterSet.Config as cms

diLeptonFilter = cms.EDFilter("DiLeptonFilter", 
  primarySrc = cms.InputTag("cleanLayer1Muons"),
  secondarySrc = cms.InputTag("cleanLayer1Muons"),
  sameSign = cms.bool(True),
  matching = cms.bool(False),
  method = cms.string('inclusive'),
  minDR = cms.double(0.0000),
  minDpt = cms.double(0.00)
)
