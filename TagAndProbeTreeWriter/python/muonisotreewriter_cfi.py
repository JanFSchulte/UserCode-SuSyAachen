import FWCore.ParameterSet.Config as cms

muonisotreewriter = cms.EDAnalyzer('MuonIsoTreeWriter',
    src = cms.InputTag("cleanLayer1Muons"),
)
