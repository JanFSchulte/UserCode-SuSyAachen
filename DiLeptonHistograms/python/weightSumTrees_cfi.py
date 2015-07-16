import FWCore.ParameterSet.Config as cms


weightSumTrees = cms.EDAnalyzer("weightSumTrees",
   genInfo = cms.InputTag("generator"),		   
)

