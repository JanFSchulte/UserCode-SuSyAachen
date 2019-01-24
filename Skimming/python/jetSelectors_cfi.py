import FWCore.ParameterSet.Config as cms

patJetFlagFilter = cms.EDFilter("PATJetFlagSelector", filter = cms.bool(False),
  src = cms.InputTag("basicJets"),
  cut = cms.string("")
)

patPFJetIDFilter = cms.EDFilter("PATPFJetIDSelector", 
  filter = cms.bool(True),
  src = cms.InputTag("basicJets"),
  version = cms.string('FIRSTDATA'),
  quality = cms.string('LOOSE')
)

jetMuonCleaner = cms.EDProducer('jetMuonCleaner',
    src = cms.InputTag("basicJets"),
    leptSrc = cms.InputTag("basicMuons"),
    dRJetLepton = cms.double(0.4)
)

jetElectronCleaner = cms.EDProducer('jetElectronCleaner',
    src = cms.InputTag("basicJets"),
    leptSrc = cms.InputTag("basicElectron"),
    dRJetLepton = cms.double(0.4)
)


