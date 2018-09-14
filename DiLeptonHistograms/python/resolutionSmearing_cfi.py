import FWCore.ParameterSet.Config as cms

resolutionSmearing = cms.PSet(   
    ResolutionSmearingFileName = cms.string('resolutionHistos.root'),    
    DoSmearing = cms.bool(True)
)
