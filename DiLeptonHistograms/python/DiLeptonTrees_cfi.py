import FWCore.ParameterSet.Config as cms

DiLeptonTrees = cms.EDAnalyzer("DiLeptonTrees",
   electrons = cms.InputTag("triggerMatchedPatElectronsPF"),
   muons = cms.InputTag("triggerMatchedPatMuonsPF"),
   taus = cms.InputTag("triggerMatchedPatTausPF"),
   jets = cms.InputTag("triggerMatchedPatJetsPF"),
   met = cms.InputTag("patMETsPF"),
   susyVars = cms.VPSet()

)

DiLeptonTreesmSugra = cms.EDAnalyzer("DiLeptonTrees",
   electrons = cms.InputTag("triggerMatchedPatElectronsPF"),
   muons = cms.InputTag("triggerMatchedPatMuonsPF"),
   taus = cms.InputTag("triggerMatchedPatTausPF"),
   jets = cms.InputTag("triggerMatchedPatJetsPF"),
   met = cms.InputTag("patMETsPF"),
   susyVars = cms.VPSet(
       cms.PSet(var = cms.string("susyScanA0"), type = cms.string("float")),
       cms.PSet(var = cms.string("susyScanM0"), type = cms.string("float")),
       cms.PSet(var = cms.string("susyScanM12"), type = cms.string("float")),
       cms.PSet(var = cms.string("susyScantanbeta"), type = cms.string("float")),
       cms.PSet(var = cms.string("susyScanCrossSection"), type = cms.string("float")),
       cms.PSet(var = cms.string("susyScanRun"), type = cms.string("float")),
       cms.PSet(var = cms.string("susyScanMu"), type = cms.string("int"))
       )

)
