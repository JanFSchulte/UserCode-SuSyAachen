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
           beamSpotSource  = cms.InputTag("offlinePrimaryVertices") #offlinePrimeryVertices
)


PATElectronIDSelector = cms.EDFilter("PATElectronIDSelector", 
  filter = cms.bool(True),
  src = cms.InputTag("cleanLayer1Electrons"),
  idMapSource  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V1-miniAOD-standalone-loose") #offlinePrimeryVertices
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

from SuSyAachen.Skimming.electronSelection_cff import noMatchedConversionsElectrons
PATElectronNoMatchedConversionSelector = noMatchedConversionsElectrons.clone(
    )


patMatchedElectronSelector = cms.EDFilter("PATElectronMatchedSelector", 
   src = cms.InputTag("cleanLayer1Electrons"),
   pdgId = cms.int32(11),
   status = cms.uint32(1),
   autoCharge = cms.bool(True)
)

from SuSyAachen.Skimming.electronSelection_cff import pfElectronProducer
pfPatElectrons = pfElectronProducer.clone(
    src = cms.InputTag("cleanLayer1Electrons"),
)

from SuSyAachen.Skimming.electronSelection_cff import bJetElectronProducer
bJetElectrons = bJetElectronProducer.clone(
    src = cms.InputTag("cleanLayer1Electrons"),
    jetSrc = cms.InputTag("cleanLayer1Jets")
)

from SuSyAachen.Skimming.electronSelection_cff import electronMuonCleaner
muonCleanElectrons = electronMuonCleaner.clone(
    src = cms.InputTag("cleanLayer1Electrons")
)



from SuSyAachen.Skimming.countFilter_cfi import candViewCountFilter
countSelector = candViewCountFilter.clone(
    src = cms.InputTag("cleanLayer1Electrons")
)
