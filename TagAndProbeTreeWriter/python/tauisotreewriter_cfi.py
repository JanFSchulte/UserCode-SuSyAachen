import FWCore.ParameterSet.Config as cms
from SuSyAachen.DiLeptonHistograms.vertexWeights_cfi import vertexWeights as vertexWeightPars
from SuSyAachen.TagAndProbeTreeWriter.isolationFunctor_cfi import isolationDefinitions

tauisotreewriter = cms.EDAnalyzer('TauIsoTreeWriter',
src = cms.InputTag("TaNCTaus"),
jets = cms.InputTag("triggerMatchedPatJetsPF"),
met = cms.InputTag("patMETsPF"),
vertices = cms.InputTag("offlinePrimaryVertices"),

useTauExtensions = cms.bool(True),
useSecondLeptonExtensions = cms.bool(True),
useMcInfo = cms.bool(True),
vertexWeights = vertexWeightPars,
isolationDefinitions = isolationDefinitions,

)

tauisotreewriterWithSecondLepton = tauisotreewriter.clone(
   secondLeptonElectronSrc = cms.InputTag("selectedPatElectrons"),
   secondLeptonMuonSrc = cms.InputTag("selectedPatMuons"),
   secondLeptonTauSrc = cms.InputTag("selectedPatTaus"),
   secondLeptonMinDeltaR = cms.double(0.15)
)
