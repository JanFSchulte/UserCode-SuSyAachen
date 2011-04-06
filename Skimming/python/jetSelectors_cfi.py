import FWCore.ParameterSet.Config as cms

patJetCountFilter = cms.EDFilter("PATJetCountFilter", filter = cms.bool(False),
  src = cms.InputTag("basicJets"),
  minNumber = cms.uint32(1),
  maxNumber = cms.uint32(99999999),
  cut = cms.string(""),
)


candViewCountFilter = cms.EDFilter("CandViewJetCountFilter", filter = cms.bool(False),
  src = cms.InputTag("basicJets"),
  minNumber = cms.uint32(1),
  maxNumber = cms.uint32(99999999),
  cut = cms.string(""),
)

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

patJetIDFilter = cms.EDFilter("PATJetIDSelector", 
  filter = cms.bool(True),
  src = cms.InputTag("basicJets"),
  version = cms.string('PURE09'),
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

resCorrectedJetProducer = cms.EDProducer('ResCorrJetsProducer',
  src = cms.InputTag("basicJets"),
  jetCorrections = cms.string("ak5PFJetsL2L3Residual")
)

unCorrectedJetProducer = cms.EDProducer('UnCorrJetsProducer',
  src = cms.InputTag("basicJets"),
)
