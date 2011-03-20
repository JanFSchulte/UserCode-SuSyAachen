import FWCore.ParameterSet.Config as cms

electronisotreewriter = cms.EDAnalyzer('ElectronIsoTreeWriter',
src = cms.InputTag("cleanLayer1Electrons"),
jets = cms.InputTag("triggerMatchedPatJetsPF"),
met = cms.InputTag("patMETsPF")
)

electronisotreewriterWithSecondLepton = electronisotreewriter.clone(
   secondLeptonElectronSrc = cms.InputTag("selectedPatElectrons"),
   secondLeptonMuonSrc = cms.InputTag("selectedPatMuons"),
   secondLeptonTauSrc = cms.InputTag("selectedPatTaus"),
   secondLeptonMinDeltaR = cms.double(0.01)
)
