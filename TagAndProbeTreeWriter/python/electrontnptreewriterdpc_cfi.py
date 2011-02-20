import FWCore.ParameterSet.Config as cms

electrontnptreewriterdpc = cms.EDAnalyzer('ElectronTnPTreeWriterDPC',
mcInfoAvailable = cms.untracked.bool(False),
cut_TnPlowInvM = cms.untracked.double(50.),
cut_TnPhighInvM = cms.untracked.double(120.),

mcSource = cms.InputTag("genParticles"),
tnpPairsSource = cms.InputTag("EE"),
passProbeSource = cms.InputTag("cleanLayer1Electrons"),
jetSource = cms.InputTag("cleanLayer1JetsAK5"),
vertexSource = cms.InputTag("offlinePrimaryVertices")
)

electrontracktnptreewriterdpc = cms.EDAnalyzer('ElectronTrackTnPTreeWriterDPC',
mcInfoAvailable = cms.untracked.bool(False),
cut_TnPlowInvM = cms.untracked.double(50.),
cut_TnPhighInvM = cms.untracked.double(120.),

mcSource = cms.InputTag("genParticles"),
tnpPairsSource = cms.InputTag("EE"),
passProbeSource = cms.InputTag("cleanLayer1Electrons"),
jetSource = cms.InputTag("cleanLayer1JetsAK5"),
vertexSource = cms.InputTag("offlinePrimaryVertices")
)

electronisotnptreewriterdpc = cms.EDAnalyzer('ElectronIsoTnPTreeWriterDPC',
mcInfoAvailable = cms.untracked.bool(False),
cut_TnPlowInvM = cms.untracked.double(50.),
cut_TnPhighInvM = cms.untracked.double(120.),

mcSource = cms.InputTag("genParticles"),
tnpPairsSource = cms.InputTag("EE"),
passProbeSource = cms.InputTag("cleanLayer1Electrons"),
jetSource = cms.InputTag("cleanLayer1JetsAK5"),
vertexSource = cms.InputTag("offlinePrimaryVertices")
)
