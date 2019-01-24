import FWCore.ParameterSet.Config as cms


metFilterNamesData2017=cms.untracked.vstring(
    "Flag_goodVertices",
    "Flag_globalSuperTightHalo2016Filter",
    "Flag_HBHENoiseFilter",
    "Flag_HBHENoiseIsoFilter",
    "Flag_EcalDeadCellTriggerPrimitiveFilter",
    "Flag_BadPFMuonFilter",
    "Flag_BadChargedCandidateFilter",
    "Flag_eeBadScFilter",
    "Flag_ecalBadCalibFilter"
)

metFilterNamesFullSim2017=cms.untracked.vstring(
    "Flag_goodVertices",
    "Flag_globalSuperTightHalo2016Filter",
    "Flag_HBHENoiseFilter",
    "Flag_HBHENoiseIsoFilter",
    "Flag_EcalDeadCellTriggerPrimitiveFilter",
    "Flag_BadPFMuonFilter",
    "Flag_BadChargedCandidateFilter",
    "Flag_ecalBadCalibFilter"
)

metFilterNamesFastSim2017=cms.untracked.vstring(                                      
    "Flag_goodVertices",              
    "Flag_HBHENoiseFilter",
    "Flag_HBHENoiseIsoFilter",
    "Flag_EcalDeadCellTriggerPrimitiveFilter",
    "Flag_BadPFMuonFilter",
    "Flag_BadChargedCandidateFilter",
    "Flag_ecalBadCalibFilter"                      
)
