import FWCore.ParameterSet.Config as cms

defaultSelector = cms.EDFilter("PATElectronSelector", 
     filter = cms.bool(True),
           src = cms.InputTag("cleanLayer1Electrons"),
           cut = cms.string('abs( eta ) <= 2.5 & pt >= 10')#GeV
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
           beamSpotSource  = cms.InputTag("offlinePrimaryVertices") #offlinePrimaryVertices
           
)


PATElectronMVAIDSelector = cms.EDFilter("PATElectronMVAIDSelector", 
  filter = cms.bool(True),
  workingPointCentralBarrelHighPt = cms.double(0.52),
  workingPointCentralBarrelLowPt = cms.double(0.2),
  workingPointCentralBarrelLowPtLinear = cms.double(0.032),
  workingPointOuterBarrelHighPt = cms.double(0.11),
  workingPointOuterBarrelLowPt = cms.double(0.1),
  workingPointOuterBarrelLowPtLinear = cms.double(0.025),
  workingPointEndcapHighPt = cms.double(0.32),
  workingPointEndcapLowPt = cms.double(-0.1),
  workingPointEndcapLowPtLinear = cms.double(0.028),
  src = cms.InputTag("cleanLayer1Electrons"),
  #idMapSource  = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV1Values") #offlinePrimeryVertices
  idMapSource  = cms.string("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV1Values") #offlinePrimeryVertices
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
