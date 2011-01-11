import FWCore.ParameterSet.Config as cms

tauisotreewriter = cms.EDAnalyzer('TauIsoTreeWriter',
src = cms.InputTag("TaNCTaus"),
jets = cms.InputTag("triggerMatchedPatJetsPF"),
met = cms.InputTag("patMETsPF")
)
