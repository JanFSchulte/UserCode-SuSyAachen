import FWCore.ParameterSet.Config as cms

BTagCalibration = cms.PSet(
    CSVFullSimTagger = cms.string('CSVv2'),
    CSVFastSimTagger = cms.string('CSV'),
    CSVFullSimFileName = cms.string('CSVv2.csv'),    
    CSVFastSimFileName = cms.string('CSV_13TEV_Combined_20_11_2015.csv'),    
    #~ CSVFullSimFileName = cms.string('/afs/cern.ch/user/c/cschomak/public/CSVv2.csv'),    
    #~ CSVFastSimFileName = cms.string('/afs/cern.ch/user/c/cschomak/public/CSV_13TEV_Combined_20_11_2015.csv'),    
)
