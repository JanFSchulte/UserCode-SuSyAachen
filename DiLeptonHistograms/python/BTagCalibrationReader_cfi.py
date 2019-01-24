import FWCore.ParameterSet.Config as cms

BTagCalibrationReaderPars2017 = cms.PSet(
    measurementType_bJets = cms.string('comb'),
    measurementType_cJets = cms.string('comb'),
    measurementType_lightJets = cms.string('incl'),
    measurementType_fastSim = cms.string('fastsim'),
    
)


BTagCalibrationReaderPars2016 = cms.PSet(
    measurementType_bJets = cms.string('comb'),
    measurementType_cJets = cms.string('comb'),
    measurementType_lightJets = cms.string('incl'),
    measurementType_fastSim = cms.string('fastsim'),
    
)
