import FWCore.ParameterSet.Config as cms

BTagCalibrationReader = cms.PSet(
    measurementType_bJets = cms.string('mujets'),
    measurementType_cJets = cms.string('mujets'),
    measurementType_lightJets = cms.string('comb'),
    measurementType_fastSim = cms.string('fastsim'),
    
)
