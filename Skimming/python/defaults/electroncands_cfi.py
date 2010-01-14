import FWCore.ParameterSet.Config as cms

defaultSelector = cms.EDFilter("CandSelector",
        src = cms.InputTag("allSuperClusters"),
        cut = cms.string('et  > 20.0 & ((abs( eta ) < 1.4442) | (abs( eta ) > 1.560 & abs( eta ) < 2.5) )')
)

