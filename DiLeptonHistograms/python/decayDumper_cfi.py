import FWCore.ParameterSet.Config as cms

decayPrinter = cms.EDAnalyzer("ParticleDecayDrawer",
   src=cms.InputTag("genParticles"),
   printP4=cms.untracked.bool(False),
   printPtEtaPhi=cms.untracked.bool(False),
   printVertex=cms.untracked.bool(False)
)

treePrinter = cms.EDAnalyzer("ParticleTreeDrawer",
    src=cms.InputTag("genParticles"),
    printP4=cms.untracked.bool(False),
    printPtEtaPhi=cms.untracked.bool(False),
    printVertex=cms.untracked.bool(False),
    printStatus=cms.untracked.bool(False),
    printIndex=cms.untracked.bool(False)
)
