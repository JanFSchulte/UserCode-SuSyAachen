import FWCore.ParameterSet.Config as cms

leptonCounter  = cms.EDAnalyzer('LeptonCounter', 
  subDir = cms.string('leptonCounter'),
  method = cms.string('exclusive'),
  electronSource = cms.InputTag("cleanLayer1Electrons"),
  muonSource = cms.InputTag("cleanLayer1Muons"),
  tauSource = cms.InputTag("cleanLayer1Taus"),
)
