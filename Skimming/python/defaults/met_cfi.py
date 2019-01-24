import FWCore.ParameterSet.Config as cms

defaultSelector = cms.EDFilter("PATMETSelector",
           filter=cms.bool(True),
           src=cms.InputTag("layer1METsAK5"),
           cut=cms.string("pt > 100")
)

