import FWCore.ParameterSet.Config as cms

defaultSelector = cms.EDFilter("GenericTrackSelector", 
	   filter = cms.bool(True),
           src = cms.InputTag("generalTracks"),
           cut = cms.string('abs( eta ) <= 2.0 & pt >= 10')#GeV
)

