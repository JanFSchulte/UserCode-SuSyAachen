import FWCore.ParameterSet.Config as cms

HLTTreeWriter = cms.EDAnalyzer('HLTTreeWriter',
onlineJets = cms.InputTag("hltJetForHT"),
onlineJetsMHT = cms.InputTag("hltJetForMHT"),
jets = cms.InputTag("triggerMatchedPatJetsPF"),
met = cms.InputTag("patMETsPF"),
jetsNoEta = cms.InputTag("triggerMatchedPatJetsPF"),
)
