import FWCore.ParameterSet.Config as cms

vertexWeights = cms.PSet(
    mcFile = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/PU_dist.root'),
    mcName = cms.string('pileup'),
    dataFile = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/PU_dist.root'),
    dataName = cms.string('pileup'),
)
