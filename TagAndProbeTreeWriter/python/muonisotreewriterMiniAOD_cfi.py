import FWCore.ParameterSet.Config as cms
from SuSyAachen.DiLeptonHistograms.vertexWeights_cfi import vertexWeights as vertexWeightPars
from SuSyAachen.TagAndProbeTreeWriter.isolationFunctor_cfi import isolationDefinitions

muonisotreewriterMiniAOD = cms.EDAnalyzer('MuonIsoTreeWriterMiniAOD',
    src = cms.InputTag("cleanLayer1Muons"),
    jets = cms.InputTag("triggerMatchedPatJetsPF"),
    met = cms.InputTag("patMETsPF"),
    vertices = cms.InputTag("offlinePrimaryVertices"),

    useTauExtensions = cms.bool(False),
    useMcInfo = cms.bool(True),
    vertexWeights = vertexWeightPars,
    isolationDefinitions = isolationDefinitions,
)

muonisotreewriterMiniAODWithSecondLepton = muonisotreewriterMiniAOD.clone(
   secondLeptonElectronSrc = cms.InputTag("selectedPatElectrons"),
   secondLeptonMuonSrc = cms.InputTag("selectedPatMuons"),
   secondLeptonTauSrc = cms.InputTag("selectedPatTaus"),
   secondLeptonMinDeltaR = cms.double(0.01)
)
