import FWCore.ParameterSet.Config as cms
from SuSyAachen.DiLeptonHistograms.vertexWeights_cfi import vertexWeights as vertexWeightPars
from SuSyAachen.TagAndProbeTreeWriter.isolationFunctor_cfi import isolationDefinitions

electronisotreewriterMiniAOD = cms.EDAnalyzer('ElectronIsoTreeWriterMiniAOD',
                                       
src = cms.InputTag("cleanLayer1Electrons"),
jets = cms.InputTag("triggerMatchedPatJetsPF"),
met = cms.InputTag("patMETsPF"),
vertices = cms.InputTag("offlinePrimaryVertices"),

useTauExtensions = cms.bool(False),
useMcInfo = cms.bool(True),
vertexWeights = vertexWeightPars,
isolationDefinitions = isolationDefinitions, 
)

electronisotreewriterMiniAODWithSecondLepton = electronisotreewriterMiniAOD.clone(
   secondLeptonElectronSrc = cms.InputTag("selectedPatElectrons"),
   secondLeptonMuonSrc = cms.InputTag("selectedPatMuons"),
   secondLeptonTauSrc = cms.InputTag("selectedPatTaus"),
   secondLeptonMinDeltaR = cms.double(0.01)
)
