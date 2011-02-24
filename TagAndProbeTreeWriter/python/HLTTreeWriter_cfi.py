import FWCore.ParameterSet.Config as cms

HLTTreeWriter = cms.EDAnalyzer('HLTTreeWriter',
onlineJets = cms.InputTag("hltJet30Ht"),
jets = cms.InputTag("triggerMatchedPatJetsPF"),
met = cms.InputTag("patMETsPF")

)
