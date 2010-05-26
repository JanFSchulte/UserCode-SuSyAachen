import FWCore.ParameterSet.Config as cms

defaultSelector = cms.EDFilter("VertexSelector",
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && ndof > 4 && abs(z) <= 15 && position.Rho <= 2"), 
    filter = cms.bool(True),   
)


patDefaultSelector = cms.EDFilter("PATPrimaryVertexCleaner",
                               filter = cms.bool(True),
                               src = cms.InputTag("cleanLayer1Muons"),
                               minMultiplicity = cms.uint32(0),
                               minPtSum = cms.double(0.0),
                               maxTrackEta = cms.double(999999.0),
                               maxNormChi2 = cms.double(999999.0),
                               maxDeltaR = cms.double(999999.0),
                               maxDeltaZ = cms.double(999999.0),
                               )
