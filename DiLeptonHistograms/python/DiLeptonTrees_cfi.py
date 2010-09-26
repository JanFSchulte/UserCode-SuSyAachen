import FWCore.ParameterSet.Config as cms

DiLeptonTrees = cms.EDAnalyzer("DiLeptonTrees",
   electrons = cms.InputTag("triggerMatchedPatElectronsPF"),
   muons = cms.InputTag("triggerMatchedPatMuonsPF"),
   taus = cms.InputTag("triggerMatchedPatTausPF"),

   electronMass = cms.double(0.000510998910), #GeV
   muonMass = cms.double(0.105658367), #GeV
   tauMass = cms.double(1.77682), #GeV
)
