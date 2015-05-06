import FWCore.ParameterSet.Config as cms

defaultTriggerDefinition =  cms.PSet(
    bits = cms.InputTag("TriggerResults","","HLT"),
    prescales = cms.InputTag("patTrigger"),
    objects = cms.InputTag("selectedPatTrigger"),
)
