import FWCore.ParameterSet.Config as cms

BTagCalibrationReader = cms.PSet(
    measurementType_bJets = cms.string('comb'),
    measurementType_cJets = cms.string('comb'),
    measurementType_lightJets = cms.string('incl'),
    measurementType_fastSim = cms.string('fastsim'),
    
)
