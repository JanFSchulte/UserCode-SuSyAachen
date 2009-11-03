import FWCore.ParameterSet.Config as cms

patMatchedElectronSelector = cms.EDFilter("PATElectronMatchedSelector", 
   src = cms.InputTag("cleanLayer1Electrons"),
   pdgId = cms.int32(11),
   status = cms.uint32(1),
   autoCharge = cms.bool(True)
)

patMatchedMuonSelector = cms.EDFilter("PATMuonMatchedSelector", 
   src = cms.InputTag("cleanLayer1Muons"),
   pdgId = cms.int32(13),
   status = cms.uint32(1),
   autoCharge = cms.bool(True)
)

patMatchedTauSelector = cms.EDFilter("PATTauMatchedSelector", 
   src = cms.InputTag("cleanLayer1Taus"),
   pdgId = cms.int32(15),
   status = cms.uint32(0),
   autoCharge = cms.bool(True)
)
