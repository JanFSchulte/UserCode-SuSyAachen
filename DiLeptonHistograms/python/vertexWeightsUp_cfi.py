import FWCore.ParameterSet.Config as cms

vertexWeightsUp = cms.PSet(
    mcFile = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/mc_PU_dist_Spring2016.root'),
    mcName = cms.string('mix_2016_25ns_SpringMC_PUScenarioV1_PoissonOOTPU'),
    #~ mcName = cms.string('pileup'),
    dataFile = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/PU_dist_2016_Up.root'),
    dataName = cms.string('pileup'),
    mc3DFile = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/mc_PU_3D_dist.root'),
    mc3DName = cms.string('pileup'),
    data3DFile = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/PU_3D_dist.root'),
    data3DName = cms.string('pileup'),
    doWeight = cms.bool(True),
    doWeight3D = cms.bool(False),
    fractionRunA = cms.double(0.46502),
    fractionRunB = cms.double(0.53498),
)
