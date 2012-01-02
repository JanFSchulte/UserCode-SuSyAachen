import FWCore.ParameterSet.Config as cms

defaultPdgIdDefinition =  cms.PSet(
    deltaR = cms.double(0.5),
    genSrc = cms.InputTag("genParticles")
)
