import FWCore.ParameterSet.Config as cms

defaultSelector = cms.EDFilter("HBHENoiseFilterMiniAOD",
	src = cms.InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResult"),
	applyfilter = cms.untracked.bool(True),
	debugOn = cms.untracked.bool(False),

)
