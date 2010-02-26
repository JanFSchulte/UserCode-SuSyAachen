import FWCore.ParameterSet.Config as cms

defaultSelector = cms.EDFilter("GenParticleSelector", 
        filter = cms.bool(True),
        src = cms.InputTag("genParticles"),
        cut = cms.string('abs( eta ) <= 5.0 & pt >= 5')#GeV
)
