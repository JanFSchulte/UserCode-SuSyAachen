import FWCore.ParameterSet.Config as cms

defaultSelector = cms.EDFilter("PATPrimaryVertexCleaner",
                               filter = cms.bool(True),
                               src = cms.InputTag("cleanLayer1Muons"),
                               minMultiplicity = cms.uint32(0),
                               minPtSum = cms.double(0.0),
                               maxTrackEta = cms.double(999999.0),
                               maxNormChi2 = cms.double(999999.0),
                               maxDeltaR = cms.double(999999.0),
                               maxDeltaZ = cms.double(999999.0),
                               )
