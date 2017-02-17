import FWCore.ParameterSet.Config as cms

vertexWeightsDown = cms.PSet(
    mcFile = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/mc_PU_dist_Spring2016.root'),
    mcName = cms.string('mix_2016_25ns_SpringMC_PUScenarioV1_PoissonOOTPU'),
    dataFile = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/PU_dist_2016_Down.root'),
    dataName = cms.string('pileup'),
    doWeight = cms.bool(True),
)
