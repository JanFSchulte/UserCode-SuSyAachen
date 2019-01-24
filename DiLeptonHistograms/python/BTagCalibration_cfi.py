import FWCore.ParameterSet.Config as cms

BTagCalibrationPars2017 = cms.PSet(
    CSVFullSimTagger = cms.string('CSVv2'),
    CSVFastSimTagger = cms.string('CSVv2'),
    CSVFullSimFileName = cms.string('CSVv2_94XSF_V2_B_F.csv'),    
    CSVFastSimFileName = cms.string('fastsim_csvv2_ttbar_26_1_2017.csv'),    
)
