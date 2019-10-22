import FWCore.ParameterSet.Config as cms


BTagCalibrationPars2016 = cms.PSet(
    CSVFullSimTagger = cms.string('DeepCSV'),
    CSVFastSimTagger = cms.string('DeepCSV'),
    CSVFullSimFileName = cms.string('DeepCSV_2016LegacySF_V1.csv'),    
    CSVFastSimFileName = cms.string('deepcsv_13TEV_16SL_18_3_2019.csv'),    
)

BTagCalibrationPars2017 = cms.PSet(
    CSVFullSimTagger = cms.string('DeepCSV'),
    CSVFastSimTagger = cms.string('DeepCSV'),
    CSVFullSimFileName = cms.string('DeepCSV_94XSF_V3_B_F.csv'),    
    CSVFastSimFileName = cms.string('deepcsv_13TEV_17SL_18_3_2019.csv'),    
)

BTagCalibrationPars2018 = cms.PSet(
    CSVFullSimTagger = cms.string('DeepCSV'),
    CSVFastSimTagger = cms.string('DeepCSV'),
    CSVFullSimFileName = cms.string('DeepCSV_102XSF_V1.csv'),    
    CSVFastSimFileName = cms.string('deepcsv_13TEV_18SL_7_5_2019.csv'),    
)
