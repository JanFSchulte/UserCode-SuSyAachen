import FWCore.ParameterSet.Config as cms
from SuSyAachen.DiLeptonHistograms.vertexWeights_cfi import vertexWeights as vertexWeightPars

tauisotreewriter = cms.EDAnalyzer('TauIsoTreeWriter',
src = cms.InputTag("TaNCTaus"),
jets = cms.InputTag("triggerMatchedPatJetsPF"),
met = cms.InputTag("patMETsPF"),
vertices = cms.InputTag("offlinePrimaryVertices"),
vertexWeights = vertexWeightPars,
)

tauisotreewriterWithSecondLepton = tauisotreewriter.clone(
   secondLeptonElectronSrc = cms.InputTag("selectedPatElectrons"),
   secondLeptonMuonSrc = cms.InputTag("selectedPatMuons"),
   secondLeptonTauSrc = cms.InputTag("selectedPatTaus"),
   secondLeptonMinDeltaR = cms.double(0.15)
)
