import FWCore.ParameterSet.Config as cms

defaultSelector = cms.EDFilter("InverseMETFilter",
	src = cms.InputTag("beamHaloFilter"),
	applyfilter = cms.untracked.bool(True),
	debugOn = cms.untracked.bool(False),

)
