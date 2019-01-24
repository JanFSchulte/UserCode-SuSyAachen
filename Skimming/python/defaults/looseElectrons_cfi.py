import FWCore.ParameterSet.Config as cms

defaultSelector = cms.EDFilter("PATElectronSelector", 
     filter = cms.bool(True),
           src = cms.InputTag("cleanLayer1Electrons"),
           cut = cms.string('abs( eta ) <= 2.0 & pt >= 10')#GeV
)

PATElectronD0BSSelector = cms.EDFilter("PATElectronD0BSSelector", 
           filter = cms.bool(True),
           src = cms.InputTag("cleanLayer1Electrons"),
           d0Min = cms.double(-999999.),
           d0Max = cms.double(0.2),
           dZMin = cms.double(-999999.),
           dZMax = cms.double(1.),
           beamSpotSource  = cms.InputTag("offlineBeamSpot") #offlinePrimeryVertices
)

PATElectronD0PVSelector = cms.EDFilter("PATElectronD0PVSelector", 
           filter = cms.bool(True),
           src = cms.InputTag("cleanLayer1Electrons"),
           d0Min = cms.double(-999999.),
           d0MaxEB = cms.double(0.2),
           d0MaxEE = cms.double(0.2),
           dZMin = cms.double(-999999.),
           dZMaxEB = cms.double(1.),
           dZMaxEE = cms.double(1.),
           SIP3DMin = cms.double(-99999.),
           SIP3DMax = cms.double(8.),
           beamSpotSource  = cms.InputTag("offlinePrimaryVertices") #offlinePrimeryVertices
           
)


PATElectronMVAIDSelector = cms.EDFilter("PATElectronMVAIDSelector", 
  filter = cms.bool(True),
  workingPointCentralBarrelHighPt = cms.double(-0.96),
  workingPointCentralBarrelLowPt = cms.double(-0.86),
  workingPointCentralBarrelLowPtLinear = cms.double(-0.86),
  workingPointOuterBarrelHighPt = cms.double(-0.96),
  workingPointOuterBarrelLowPt = cms.double(-0.85),
  workingPointOuterBarrelLowPtLinear = cms.double(-0.85),
  workingPointEndcapHighPt = cms.double(-0.95),
  workingPointEndcapLowPt = cms.double(-0.81),
  workingPointEndcapLowPtLinear = cms.double(-0.81),
  src = cms.InputTag("cleanLayer1Electrons"),
  idMapSource  = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV1Values") #offlinePrimeryVertices
)

PATElectronConversionSelector = cms.EDFilter("PATElectronConversionSelector",
           filter = cms.bool(True),
           src = cms.InputTag("cleanLayer1Electrons"),
           convInfo = cms.bool(False),
           maxDistance = cms.double(0.02),
           maxCotangentTheta = cms.double(0.02),
           minLostHits = cms.int32(0),
           maxLostHits = cms.int32(1)
)

from SuSyAachen.Skimming.electronSelection_cff import effectiveAreaIsoElectrons
PATElectronEffectiveAreaSelector = effectiveAreaIsoElectrons.clone(
)
from SuSyAachen.Skimming.electronSelection_cff import isoElectrons
PATElectronIsolationSelector = isoElectrons.clone(
)

from SuSyAachen.Skimming.electronSelection_cff import noMatchedConversionsElectrons
PATElectronNoMatchedConversionSelector = noMatchedConversionsElectrons.clone(
    )



