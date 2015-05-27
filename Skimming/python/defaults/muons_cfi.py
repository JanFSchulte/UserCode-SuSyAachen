import FWCore.ParameterSet.Config as cms

defaultSelector = cms.EDFilter("PATMuonSelector", 
	   filter = cms.bool(True),
           src = cms.InputTag("cleanLayer1Muons"),
           cut = cms.string('abs( eta ) <= 2.0 & pt >= 10')#GeV
)

PATMuonD0BSSelector = cms.EDFilter("PATMuonD0BSSelector", 
           filter = cms.bool(True),
           src = cms.InputTag("cleanLayer1Muons"),
           d0Min = cms.double(-999999.),
           d0Max = cms.double(0.2),
           dZMin = cms.double(-999999.),
           dZMax = cms.double(1.),
           beamSpotSource  = cms.InputTag("offlineBeamSpot") #offlinePrimeryVertices
)

PATMuonD0PVSelector = cms.EDFilter("PATMuonD0PVSelector", 
           filter = cms.bool(True),
           src = cms.InputTag("cleanLayer1Muons"),
           d0Min = cms.double(-999999.),
           d0Max = cms.double(0.2),
           dZMin = cms.double(-999999.),
           dZMax = cms.double(1.),
           beamSpotSource  = cms.InputTag("offlinePrimaryVertices") #offlinePrimeryVertices
)

PATMuonTightIDSelector = cms.EDFilter("PATMuonTightIDSelector", 
           filter = cms.bool(True),
           src = cms.InputTag("cleanLayer1Muons"),
           vertexSource  = cms.InputTag("offlineSlimmedPrimaryVertices") #offlinePrimeryVertices
)

patMatchedMuonSelector = cms.EDFilter("PATMuonMatchedSelector", 
   src = cms.InputTag("cleanLayer1Muons"),
   pdgId = cms.int32(13),
   status = cms.uint32(1),
   autoCharge = cms.bool(True)
)

from SuSyAachen.Skimming.muonSelection_cff import bJetMuonProducer
bJetMuons = bJetMuonProducer.clone(
    src = cms.InputTag("cleanLayer1Muons"),
    jetSrc = cms.InputTag("cleanLayer1Jets")
)

from SuSyAachen.Skimming.countFilter_cfi import candViewCountFilter
countSelector = candViewCountFilter.clone(
    src = cms.InputTag("cleanLayer1Muons")
)
from SuSyAachen.Skimming.muonSelection_cff import isoMuons
PATMuonIsolationSelector = isoMuons.clone(
)

