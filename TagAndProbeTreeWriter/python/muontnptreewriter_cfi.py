import FWCore.ParameterSet.Config as cms

muontnptreewriter = cms.EDAnalyzer('MuonTnPTreeWriter',
mcInfoAvailable = cms.untracked.bool(False),

mcSource = cms.InputTag("genParticles"),
tagSource = cms.InputTag("cleanLayer1Muons"),
probeSource = cms.InputTag("generalTracks"),
jetSource = cms.InputTag("cleanLayer1JetsAK5")
)