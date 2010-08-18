import FWCore.ParameterSet.Config as cms

electrontnptreewriter = cms.EDAnalyzer('ElectronTnPTreeWriter',
mcInfoAvailable = cms.untracked.bool(False),
cut_TnPlowInvM = cms.untracked.double(50.),
cut_TnPhighInvM = cms.untracked.double(120.),

mcSource = cms.InputTag("genParticles"),
tagSource = cms.InputTag("cleanLayer1Electrons"),
probeSource = cms.InputTag("cleanLayer1Electrons"),
passProbeSource = cms.InputTag("cleanLayer1Electrons"),
jetSource = cms.InputTag("cleanLayer1JetsAK5")
)

electronisotnptreewriter = cms.EDAnalyzer('ElectronIsoTnPTreeWriter',
mcInfoAvailable = cms.untracked.bool(False),
cut_TnPlowInvM = cms.untracked.double(50.),
cut_TnPhighInvM = cms.untracked.double(120.),

mcSource = cms.InputTag("genParticles"),
tagSource = cms.InputTag("cleanLayer1Electrons"),
probeSource = cms.InputTag("cleanLayer1Electrons"),
passProbeSource = cms.InputTag("cleanLayer1Electrons"),
jetSource = cms.InputTag("cleanLayer1JetsAK5")
)
