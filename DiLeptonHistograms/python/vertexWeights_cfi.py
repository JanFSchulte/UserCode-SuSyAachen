import FWCore.ParameterSet.Config as cms

vertexWeightsPars2017 = cms.PSet(
    mcFile = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2017/mc_PU_dist_Fall2017.root'),
    mcName = cms.string('mix_2017_25ns_WinterMC_PUScenarioV1_PoissonOOTPU'),
    dataFile = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2017/PU_dist_2017.root'),
    dataName = cms.string('pileup'),
    doWeight = cms.bool(True),
    verbosity = cms.int32(0),
)

vertexWeightsPars2016 = cms.PSet(
    mcFile = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2016/mc_PU_dist_Summer2016.root'),
    mcName = cms.string('mix_2016_25ns_Moriond17MC_PoissonOOTPU'),
    dataFile = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2016/PU_dist_2016.root'),
    dataName = cms.string('pileup'),
    doWeight = cms.bool(True),
    verbosity = cms.int32(0),
)
