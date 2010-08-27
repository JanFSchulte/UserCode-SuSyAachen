import FWCore.ParameterSet.Config as cms

patLeptonCountFilter = cms.EDFilter("PATLeptonCountFilter",
  filter = cms.bool(False),
  muonSource = cms.InputTag("basicMuons"),
  electronSource = cms.InputTag("basicElectrons"),
  tauSource = cms.InputTag("basicTaus"),

  countElectrons = cms.bool(True),
  countMuons = cms.bool(True),
  countTaus = cms.bool(True),
                                    
  minNumber = cms.uint32(1),
  maxNumber = cms.uint32(999999),
  cut = cms.string(""),
)

bJetElectronProducer = cms.EDProducer('bJetElectronProducer',
          src = cms.InputTag("basicElectrons"),
          jetSrc = cms.InputTag("basicJets"),
          dRJetLepton = cms.double(0.2),
          dPhiOppositeJetLepton = cms.double(2.7),
          user_bJetAlgo = cms.untracked.string("trackCountingHighPurBJetTags"),
          user_bTagDiscriminator = cms.untracked.double(3.),
          )

bJetMuonProducer = cms.EDProducer('bJetMuonProducer',
          src = cms.InputTag("basicMuons"),
          jetSrc = cms.InputTag("basicJets"),
          dRJetLepton = cms.double(0.2),
          dPhiOppositeJetLepton = cms.double(2.7),
          user_bJetAlgo = cms.untracked.string("trackCountingHighPurBJetTags"),
          user_bTagDiscriminator = cms.untracked.double(3.),
          )

