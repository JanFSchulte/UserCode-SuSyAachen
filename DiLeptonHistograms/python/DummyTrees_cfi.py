import FWCore.ParameterSet.Config as cms


DummyTrees = cms.EDAnalyzer('DummyTrees',
electrons = cms.InputTag("cleanLayer1Electrons"),
muons = cms.InputTag("cleanLayer1Muons"),
taus = cms.InputTag("cleanLayer1Taus"),
jets = cms.InputTag("cleanLayer1JetsAK4"),
met = cms.InputTag("cleanMETsAK4"),

bJetAlgo = cms.untracked.string("trackCountingHighEffBJetTags"),
bTagDiscriminator = cms.untracked.double(1.7)
)
