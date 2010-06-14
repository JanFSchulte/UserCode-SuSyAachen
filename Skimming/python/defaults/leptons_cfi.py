import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.selectionLayer1.leptonCountFilter_cfi import countPatLeptons
defaultSelector = countPatLeptons.clone(
    filter = cms.bool(False),
    electronSource = cms.InputTag("cleanPatElectrons"),
    muonSource     = cms.InputTag("cleanPatMuons"),
    tauSource      = cms.InputTag("cleanPatTaus"),
    countElectrons = cms.bool(True),
    countMuons     = cms.bool(True),
    countTaus      = cms.bool(False),
    minNumber = cms.uint32(0),
    maxNumber = cms.uint32(999999),
)

