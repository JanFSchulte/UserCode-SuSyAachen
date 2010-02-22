import FWCore.ParameterSet.Config as cms

MCDiLeptonHistograms = cms.EDAnalyzer("MCDiLeptonHistograms",
   genParticles = cms.InputTag("genParticles"),
   electrons = cms.InputTag("promptElectrons"),#promptElectrons"),
   muons = cms.InputTag("promptMuons"),
   taus = cms.InputTag("promptTaus"),
   tauJets = cms.InputTag("tauGenJets"),
)
