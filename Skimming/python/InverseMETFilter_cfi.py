import FWCore.ParameterSet.Config as cms

InverseMETFilter = cms.EDFilter("InverseMETFilter", 
   hcalFilterTag = cms.InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResult"),
   beamHaloTag = cms.InputTag("beamHaloFilter"),
   ecalDeadCellTag = cms.InputTag("EcalDeadCellTriggerPrimitiveFilter"),
   eeBadScTag = cms.InputTag("eeBadScFilter"),
   hcalLaserTag = cms.InputTag("hcalLaserEventFilter"),
   trackingFailureTag = cms.InputTag("trackingFailureFilter")

)
