import FWCore.ParameterSet.Config as cms

BTagCalibration = cms.PSet(
    CSVFullSimTagger = cms.string('CSVv2'),
    CSVFastSimTagger = cms.string('CSVv2'),
    CSVFullSimFileName = cms.string('CSVv2_94XSF_V2_B_F.csv'),    
    CSVFastSimFileName = cms.string('fastsim_csvv2_ttbar_26_1_2017.csv'),    
    #~ CSVFullSimFileName = cms.string('CSVv2_ichep.csv'),    
    #~ CSVFastSimFileName = cms.string('CSV_13TEV_Combined_14_7_2016.csv'),    
    #~ CSVFullSimFileName = cms.string('/afs/cern.ch/user/c/cschomak/public/CSVv2.csv'),    
    #~ CSVFastSimFileName = cms.string('/afs/cern.ch/user/c/cschomak/public/CSV_13TEV_Combined_20_11_2015.csv'),    
)
