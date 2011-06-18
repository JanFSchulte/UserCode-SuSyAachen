import FWCore.ParameterSet.Config as cms
from SuSyAachen.DiLeptonHistograms.vertexWeights_cfi import vertexWeights as vertexWeightPars

muonisotreewriter = cms.EDAnalyzer('MuonIsoTreeWriter',
    src = cms.InputTag("cleanLayer1Muons"),
    jets = cms.InputTag("triggerMatchedPatJetsPF"),
    met = cms.InputTag("patMETsPF"),
    vertices = cms.InputTag("offlinePrimaryVertices"),

    useTauExtensions = cms.bool(False),
    useMcInfo = cms.bool(True),
    vertexWeights = vertexWeightPars,
)

muonisotreewriterWithSecondLepton = muonisotreewriter.clone(
   secondLeptonElectronSrc = cms.InputTag("selectedPatElectrons"),
   secondLeptonMuonSrc = cms.InputTag("selectedPatMuons"),
   secondLeptonTauSrc = cms.InputTag("selectedPatTaus"),
   secondLeptonMinDeltaR = cms.double(0.01)
)
