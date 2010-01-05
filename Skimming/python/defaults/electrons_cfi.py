import FWCore.ParameterSet.Config as cms

defaultSelector = cms.EDFilter("PATElectronSelector", 
	   filter = cms.bool(True),
           src = cms.InputTag("cleanLayer1Electrons"),
           cut = cms.string('abs( eta ) <= 2.0 & pt >= 10')#GeV
)

PATElectronD0Selector = cms.EDFilter("PATElectronD0Selector", 
           filter = cms.bool(True),
           src = cms.InputTag("cleanLayer1Electrons"),
           d0Min = cms.double(0.2),
           beamSpotSource  = cms.InputTag("offlineBeamSpot") #offlinePrimeryVertices
)


patMatchedElectronSelector = cms.EDFilter("PATElectronMatchedSelector", 
   src = cms.InputTag("cleanLayer1Electrons"),
   pdgId = cms.int32(11),
   status = cms.uint32(1),
   autoCharge = cms.bool(True)
)

