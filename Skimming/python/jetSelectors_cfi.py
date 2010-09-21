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


resCorrectedJetProducer = cms.EDProducer('ResCorrJetsProducer',
  src = cms.InputTag("basicJets"),
  jetCorrections = cms.string("CondFormats/JetMETObjects/data/Spring10DataV2_L2L3Residual_AK5PF.txt")
)

unCorrectedJetProducer = cms.EDProducer('UnCorrJetsProducer',
  src = cms.InputTag("basicJets"),
)
