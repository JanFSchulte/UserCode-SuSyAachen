import FWCore.ParameterSet.Config as cms

metFilterNamesData2016=cms.untracked.vstring(                                      
  "Flag_HBHENoiseFilter", 
  "Flag_HBHENoiseIsoFilter",
  "Flag_globalSuperTightHalo2016Filter",   
  "Flag_goodVertices", 
  "Flag_EcalDeadCellTriggerPrimitiveFilter",   
  "Flag_eeBadScFilter",   
  "Flag_BadPFMuonFilter",                                                                            
) 

metFilterNamesFullSim2016=cms.untracked.vstring(                                      
  "Flag_HBHENoiseFilter", 
  "Flag_HBHENoiseIsoFilter",
  "Flag_globalSuperTightHalo2016Filter",   
  "Flag_goodVertices", 
  "Flag_EcalDeadCellTriggerPrimitiveFilter",   
  "Flag_BadPFMuonFilter",                   
) 


metFilterNamesData2017=cms.untracked.vstring(
    "Flag_goodVertices",
    "Flag_globalSuperTightHalo2016Filter",
    "Flag_HBHENoiseFilter",
    "Flag_HBHENoiseIsoFilter",
    "Flag_EcalDeadCellTriggerPrimitiveFilter",
    "Flag_BadPFMuonFilter",
    #"Flag_BadChargedCandidateFilter",
    "Flag_eeBadScFilter",
    #"ecalBadCalibReducedMINIAODFilter"
)

metFilterNamesFullSim2017=cms.untracked.vstring(
    "Flag_goodVertices",
    "Flag_globalSuperTightHalo2016Filter",
    "Flag_HBHENoiseFilter",
    "Flag_HBHENoiseIsoFilter",
    "Flag_EcalDeadCellTriggerPrimitiveFilter",
    "Flag_BadPFMuonFilter",
    #"Flag_BadChargedCandidateFilter",
    #"ecalBadCalibReducedMINIAODFilter"
)

metFilterNamesFastSim2017=cms.untracked.vstring(                                      
    "Flag_goodVertices",              
    "Flag_HBHENoiseFilter",
    "Flag_HBHENoiseIsoFilter",
    "Flag_EcalDeadCellTriggerPrimitiveFilter",
    "Flag_BadPFMuonFilter",
    #"Flag_BadChargedCandidateFilter",
    #"ecalBadCalibReducedMINIAODFilter"                      
)

metFilterNamesData2018=cms.untracked.vstring(
    "Flag_goodVertices",
    "Flag_globalSuperTightHalo2016Filter",
    "Flag_HBHENoiseFilter",
    "Flag_HBHENoiseIsoFilter",
    "Flag_EcalDeadCellTriggerPrimitiveFilter",
    "Flag_BadPFMuonFilter",
    #"Flag_BadChargedCandidateFilter",
    "Flag_eeBadScFilter",
    #"ecalBadCalibReducedMINIAODFilter"
)

metFilterNamesFullSim2018=cms.untracked.vstring(
    "Flag_goodVertices",
    "Flag_globalSuperTightHalo2016Filter",
    "Flag_HBHENoiseFilter",
    "Flag_HBHENoiseIsoFilter",
    "Flag_EcalDeadCellTriggerPrimitiveFilter",
    "Flag_BadPFMuonFilter",
    #"Flag_BadChargedCandidateFilter",
    #"ecalBadCalibReducedMINIAODFilter"
)

metFilterNamesFastSim2018=cms.untracked.vstring(                                      
    "Flag_goodVertices",              
    "Flag_HBHENoiseFilter",
    "Flag_HBHENoiseIsoFilter",
    "Flag_EcalDeadCellTriggerPrimitiveFilter",
    "Flag_BadPFMuonFilter",
    #"Flag_BadChargedCandidateFilter",
    #"ecalBadCalibReducedMINIAODFilter"                      
)
