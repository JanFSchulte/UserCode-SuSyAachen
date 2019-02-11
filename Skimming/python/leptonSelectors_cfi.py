import FWCore.ParameterSet.Config as cms

patLeptonCountFilter = cms.EDFilter("PATLeptonCountFilter",
  filter = cms.bool(False),
  muonSource = cms.InputTag("basicMuons"),
  electronSource = cms.InputTag("basicElectrons"),
  tauSource = cms.InputTag("basicTaus"),

  countElectrons = cms.bool(True),
  countMuons = cms.bool(True),
  countTaus = cms.bool(True),
                                    
  minNumber = cms.uint32(1),
  maxNumber = cms.uint32(999999),
  cut = cms.string(""),
)
