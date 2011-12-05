import FWCore.ParameterSet.Config as cms
from SuSyAachen.DiLeptonHistograms.vertexWeights_cfi import vertexWeights as vertexWeightPars


WeightSummer = cms.EDAnalyzer("WeightSummer",
   vertexWeights = vertexWeightPars,
)
