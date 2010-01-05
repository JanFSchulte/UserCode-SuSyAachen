import FWCore.ParameterSet.Config as cms

defaultSelector = cms.EDFilter("PATMuonSelector", 
	   filter = cms.bool(True),
           src = cms.InputTag("cleanLayer1Muons"),
           cut = cms.string('abs( eta ) <= 2.0 & pt >= 10')#GeV
)

PATMuonD0Selector = cms.EDFilter("PATMuonD0Selector", 
           filter = cms.bool(True),
           src = cms.InputTag("cleanLayer1Muons"),
           d0Min = cms.double(0.2),
           beamSpotSource  = cms.InputTag("offlineBeamSpot") #offlinePrimeryVertices
)

patMatchedMuonSelector = cms.EDFilter("PATMuonMatchedSelector", 
   src = cms.InputTag("cleanLayer1Muons"),
   pdgId = cms.int32(13),
   status = cms.uint32(1),
   autoCharge = cms.bool(True)
)

