import FWCore.ParameterSet.Config as cms


BTagCalibrationPars2016 = cms.PSet(
    CSVFullSimTagger = cms.string('CSVv2'),
    CSVFastSimTagger = cms.string('CSVv2'),
    CSVFullSimFileName = cms.string('CSVv2_Moriond17_B_H.csv'),    
    CSVFastSimFileName = cms.string('fastsim_csvv2_ttbar_26_1_2017.csv'),    
)

BTagCalibrationPars2017 = cms.PSet(
    CSVFullSimTagger = cms.string('DeepCSV'),
    CSVFastSimTagger = cms.string('CSVv2'),
    CSVFullSimFileName = cms.string('DeepCSV_94XSF_V3_B_F.csv'),    
    CSVFastSimFileName = cms.string('fastsim_csvv2_ttbar_26_1_2017.csv'),    
)

BTagCalibrationPars2018 = cms.PSet(
    CSVFullSimTagger = cms.string('DeepCSV'),
    CSVFastSimTagger = cms.string('CSVv2'),
    CSVFullSimFileName = cms.string('DeepCSV_102XSF_V1.csv'),    
    CSVFastSimFileName = cms.string('fastsim_csvv2_ttbar_26_1_2017.csv'),    
)
