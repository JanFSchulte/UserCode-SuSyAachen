import FWCore.ParameterSet.Config as cms

defaultSelector = cms.EDFilter("METFilterMiniAOD",
	src = cms.InputTag("TriggerResults","PAT"),
	applyfilter = cms.untracked.bool(True),
	debugOn = cms.untracked.bool(False),

)
