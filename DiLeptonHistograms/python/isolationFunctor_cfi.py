import FWCore.ParameterSet.Config as cms
isolationDefinitions = cms.PSet(
    rhoSource = cms.InputTag("fixedGridRhoFastjetCentralNeutral"),
    candSource = cms.InputTag("packedPFCandidates"),
    ) 
