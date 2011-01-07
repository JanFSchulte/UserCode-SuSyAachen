import FWCore.ParameterSet.Config as cms

TriggerAnalysis = cms.EDAnalyzer('TriggerHistograms',
    triggerSource = cms.InputTag("TriggerResults","","HLT")
)
