import FWCore.ParameterSet.Config as cms

tauisotreewriter = cms.EDAnalyzer('TauIsoTreeWriter',
src = cms.InputTag("TaNCTaus"),
jets = cms.InputTag("triggerMatchedPatJetsPF"),
met = cms.InputTag("patMETsPF"),
)

tauisotreewriterWithSecondLepton = tauisotreewriter.clone(
   secondLeptonElectronSrc = cms.InputTag("selectedPatElectrons"),
   secondLeptonMuonSrc = cms.InputTag("selectedPatMuons"),
   secondLeptonTauSrc = cms.InputTag("selectedPatTaus"),
   secondLeptonMinDeltaR = cms.double(0.15)
)
