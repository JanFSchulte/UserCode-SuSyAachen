import FWCore.ParameterSet.Config as cms

GenPromptSelector = cms.EDFilter("GenPromptSelector",
    src = cms.InputTag("genParticles"),
    promptMotherIds = cms.vint32( 15,-15, 23,-23, 24,-24),
    bsm = cms.bool(True)
)