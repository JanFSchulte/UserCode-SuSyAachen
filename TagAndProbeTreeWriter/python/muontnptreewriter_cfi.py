import FWCore.ParameterSet.Config as cms

muontnptreewriter = cms.EDAnalyzer('MuonTnPTreeWriter',
mcInfoAvailable = cms.untracked.bool(False),
cut_TnPlowInvM = cms.untracked.double(50.),
cut_TnPhighInvM = cms.untracked.double(120.),

mcSource = cms.InputTag("genParticles"),
tagSource = cms.InputTag("cleanLayer1Muons"),
probeSource = cms.InputTag("generalTracks"),
passProbeSource = cms.InputTag("cleanLayer1Muons"),
jetSource = cms.InputTag("cleanLayer1JetsAK5"),
vertexSource = cms.InputTag("offlinePrimaryVertices")
)

muonisotnptreewriter = cms.EDAnalyzer('MuonIsoTnPTreeWriter',
mcInfoAvailable = cms.untracked.bool(False),
cut_TnPlowInvM = cms.untracked.double(50.),
cut_TnPhighInvM = cms.untracked.double(120.),

mcSource = cms.InputTag("genParticles"),
tagSource = cms.InputTag("cleanLayer1Muons"),
probeSource = cms.InputTag("generalTracks"),
passProbeSource = cms.InputTag("cleanLayer1Muons"),
jetSource = cms.InputTag("cleanLayer1JetsAK5"),
vertexSource = cms.InputTag("offlinePrimaryVertices")
)
