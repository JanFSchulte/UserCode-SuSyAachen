import FWCore.ParameterSet.Config as cms

defaultSelector = cms.EDFilter("CompositeDiLeptonFilter", 
  EESrc = cms.InputTag("EE"), EE = cms.string("reject"),
  EMuSrc = cms.InputTag("EMu"), EMu = cms.string("reject"),
  MuMuSrc = cms.InputTag("MuMu"), MuMu = cms.string("reject"),
  ETauSrc = cms.InputTag("ETau"), ETau = cms.string("reject"),
  MuTauSrc = cms.InputTag("MuTau"), MuTau = cms.string("reject"),
  TauTauSrc = cms.InputTag("TauTau"), TauTau = cms.string("reject"),

  sameSign = cms.bool(True),
  channels = cms.vstring("EE","EMu", "MuMu", "ETau", "MuTau", "TauTau")
)

exclusiveDileptonFilter = defaultSelector.clone() #just for clarity

compositeProducer = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('muons muons'),
    cut = cms.string(''),
    name = cms.string('dileptonDecay'),
    roles = cms.vstring('lepton1', 'lepton2'),
    checkCharge = cms.bool(False)
)
