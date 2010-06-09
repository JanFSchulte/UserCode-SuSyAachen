import FWCore.ParameterSet.Config as cms

defaultSelector =  cms.EDProducer("PATMETSelector", 
           filter = cms.bool(True),
           src = cms.InputTag("layer1METsAK5"),
           cut = cms.string("pt > 100")
)

genMETSelector = cms.EDProducer("CandViewSelector", 
           filter = cms.bool(True),
           src = cms.InputTag("genMetTrue"),
           cut = cms.string("pt > 100")
)

from SuSyAachen.Skimming.mtFilter_cfi import mtFilter
patHtFilter = mtFilter.clone( )
