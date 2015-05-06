import FWCore.ParameterSet.Config as cms

defaultTriggerDefinition =  cms.PSet(
    triggerSrc = cms.InputTag("hltTriggerSummaryAOD","","HLT")
)
