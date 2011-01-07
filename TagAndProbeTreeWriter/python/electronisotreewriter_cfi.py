import FWCore.ParameterSet.Config as cms

electronisotreewriter = cms.EDAnalyzer('ElectronIsoTreeWriter',
src = cms.InputTag("cleanLayer1Electrons"),
jets = cms.InputTag("triggerMatchedPatJetsPF"),
met = cms.InputTag("patMETsPF")
)
