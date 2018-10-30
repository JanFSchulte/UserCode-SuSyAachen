import FWCore.ParameterSet.Config as cms

vertexWeightsDown = cms.PSet(
    mcFile = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/mc_PU_dist_Fall2017.root'),
    mcName = cms.string('mix_2017_25ns_WinterMC_PUScenarioV1_PoissonOOTPU'),
    dataFile = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/PU_dist_2017_Down.root'),
    dataName = cms.string('pileupDown'),
    doWeight = cms.bool(True),
)
