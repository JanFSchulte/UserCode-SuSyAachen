import FWCore.ParameterSet.Config as cms


pdfUncertaintyNominatorTrees = cms.EDAnalyzer("pdfUncertaintyNominatorTrees",
   pdfInfo = cms.InputTag("generator"),		   
)

