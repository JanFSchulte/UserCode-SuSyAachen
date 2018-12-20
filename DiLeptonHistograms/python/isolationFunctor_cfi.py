import FWCore.ParameterSet.Config as cms
isolationDefinitions = cms.PSet(
    rhoSource = cms.InputTag("fixedGridRhoFastjetAll"),
    candSource = cms.InputTag("packedPFCandidates"),
    ) 
