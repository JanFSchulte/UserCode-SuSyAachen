import FWCore.ParameterSet.Config as cms

GenDaughterExcluder = cms.EDFilter("GenDaughterExcluder",
    src = cms.InputTag("genParticles"),
    daughterIds = cms.vint32( 11,-11, 13,-13) #e+ e- mu+ mu-
)
