import FWCore.ParameterSet.Config as cms

diLeptonFilter = cms.EDFilter("DiLeptonFilter", 
  # use strings of size == 2; p = primary, s = secondary, t= tertiary (e.g. "pp", "ss", "ts", ...)
  combinations = cms.vstring("ps"),
  primarySrc = cms.InputTag("cleanLayer1Muons"),
  secondarySrc = cms.InputTag("cleanLayer1Muons"),
  tertiarySrc = cms.InputTag(""),# only loaded if any combination with 't' is given
  sameSign = cms.bool(True),
  matching = cms.bool(False),
  method = cms.string('inclusive'),
  minDR = cms.double(0.0000),
  minDpt = cms.double(0.00)
)
